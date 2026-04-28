#pragma once
#include "nested.h"
#include "error-heuristic.h"
#include <memory>
#include <thread>
#include <utility>
#include <vector>
#include "../relaxedPriorityQueues_874280/src/skiplist/skiplist.hpp"

// Comparador adaptado para trabajar con shared_ptr<ERegion>.
// Accede al contenido via operador-> en lugar de directamente.
struct SlComparator {
    template<typename T>
    bool operator()(const T& a, const T& b) const {
        return std::get<0>(a->extra()) < std::get<0>(b->extra());
    }
};

namespace viltrum {

template<typename Rule, typename ErrorHeuristic, typename = std::enable_if_t<is_nested<Rule>::value>>
class RegionsGeneratorAdaptiveSl {
    Rule rule;
    ErrorHeuristic error_heuristic;
    std::size_t subdivisions;

public:
    RegionsGeneratorAdaptiveSl(const Rule& r, const ErrorHeuristic& er, std::size_t subdivisions) 
        : rule(r), error_heuristic(er), subdivisions(subdivisions) {}

    template<std::size_t DIMBINS, typename F, typename Float, std::size_t DIM, typename Logger>
    auto generate(const std::array<std::size_t,DIMBINS>& bin_resolution,
		const F& f, const Range<Float,DIM>& range, Logger& logger) const {

        auto r_init = region(f, rule, range.min(), range.max());
        auto errdim_init = error_heuristic(r_init);
        using ERegion  = ExtendedRegion<decltype(r_init), decltype(errdim_init)>;
        // shared_ptr<ERegion> tiene constructor por defecto y operator= copiable,
        // por lo que satisface los requisitos internos del Skiplist (nodo NIL, extracción).
        using PERegion = std::shared_ptr<ERegion>;

        unsigned int num_threads = std::thread::hardware_concurrency();
        if (num_threads == 0) num_threads = 2;
        Skiplist<PERegion, SlComparator> sl(SlComparator(), num_threads * 8);

        logger.log_progress(std::size_t(0), subdivisions);
        sl.push(std::make_shared<ERegion>(r_init, errdim_init));

        for (std::size_t i = 0; i < subdivisions; ++i) {
            auto opt_r = sl.try_pop();
            if (!opt_r) break;

            // *opt_r es el shared_ptr; **opt_r es la ERegion apuntada.
            ERegion& current_r = **opt_r;

            auto subregions = current_r.split(f, std::get<1>(current_r.extra()));

            for (auto& sr : subregions) {
                auto errdim = error_heuristic(sr);
                sl.push(std::make_shared<ERegion>(std::move(sr), std::move(errdim)));
            }

            logger.log_progress(i, subdivisions);
        }

        // Reconvertir la Skiplist a un vector de ERegion (moviendo el contenido).
        std::vector<ERegion> final_heap;
        while (auto final_r = sl.try_pop()) {
            final_heap.push_back(std::move(**final_r));
        }

        logger.log_progress(subdivisions, subdivisions);
        return final_heap;
    }
};

template<typename Rule, typename ErrorHeuristic>
auto regions_generator_adaptive_skiplist(const Rule& r, const ErrorHeuristic& er, std::size_t subdivisions) {
    return RegionsGeneratorAdaptiveSl<Rule,ErrorHeuristic>(r,er,subdivisions);
}

} // namespace viltrum