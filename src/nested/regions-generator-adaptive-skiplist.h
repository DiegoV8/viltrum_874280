#pragma once
#include "nested.h"
#include "error-heuristic.h"
#include <memory>
#include <thread>
#include <utility>
#include <atomic>
#include <vector>
#include "../relaxedPriorityQueues_874280/src/skiplist/skiplist.hpp"

// Estructura usada para comparar dos elementos en la estructura.
struct SlComparator {
    template<typename T>
    bool operator()(const T& a, const T& b) const {
        return std::get<0>(a->extra()) > std::get<0>(b->extra());
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
    auto generate(const std::array<std::size_t, DIMBINS>& bin_resolution,
                const F& f, const Range<Float, DIM>& range, Logger& logger) const {

        auto r_init = region(f, rule, range.min(), range.max());
        auto errdim_init = error_heuristic(r_init);
        
        using ERegion = ExtendedRegion<decltype(r_init), decltype(errdim_init)>;
        using ERegionPtr = std::shared_ptr<ERegion>;

        unsigned int num_threads = std::thread::hardware_concurrency();
        if (num_threads == 0) num_threads = 2;

        Skiplist<ERegionPtr, SlComparator> sl{SlComparator(), int(num_threads)*128};
        sl.push(std::make_shared<ERegion>(r_init, errdim_init));

        std::atomic<std::size_t> completed{0};

        logger.log_progress(std::size_t(0), subdivisions);

        #pragma omp parallel
        {
            while (true) {
                std::size_t current = completed.load(std::memory_order_relaxed);
                if (current >= subdivisions) break;

                std::optional<ERegionPtr> opt_ptr = sl.try_pop();
                if (!opt_ptr) {
                    std::this_thread::yield();
                    continue;
                }

                std::size_t my_slot = completed.fetch_add(1, std::memory_order_relaxed);
                if (my_slot >= subdivisions) {
                    sl.push(*opt_ptr);
                    break;
                }

                auto subregions = (*opt_ptr)->split(f, std::get<1>((*opt_ptr)->extra()));
                for (auto& sr : subregions) {
                    sl.push(std::make_shared<ERegion>(sr, error_heuristic(sr)));
                }

                #pragma omp critical(logger_update)
                {
                    logger.log_progress(my_slot + 1, subdivisions);
                }
            }
        }

        auto final_ptrs = sl.drain();
        std::vector<ERegion> final_heap;
        final_heap.reserve(final_ptrs.size());
        for (auto& ptr : final_ptrs) {
            final_heap.push_back(std::move(*ptr));
        }

        logger.log_progress(subdivisions, subdivisions);
        return final_heap;
    }
};

template<typename Rule, typename ErrorHeuristic>
auto regions_generator_adaptive_skiplist(const Rule& r, const ErrorHeuristic& er, std::size_t subdivisions) {
    return RegionsGeneratorAdaptiveSl<Rule,ErrorHeuristic>(r,er,subdivisions);
}

}