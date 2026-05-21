#pragma once
#include "nested.h"
#include "error-heuristic.h"
#include <thread>
#include <atomic>
#include <vector>
#include "../relaxedPriorityQueues_874280/src/multiqueue/multiqueue.hpp"

struct MqComparator {
    template<typename T>
    bool operator()(const T& a, const T& b) const {
        return std::get<0>(a.extra()) < std::get<0>(b.extra());
    }
};

namespace viltrum {

template<typename Rule, typename ErrorHeuristic,
         typename = std::enable_if_t<is_nested<Rule>::value>>
class RegionsGeneratorAdaptiveMq {
    Rule rule;
    ErrorHeuristic error_heuristic;
    std::size_t subdivisions;

public:
    RegionsGeneratorAdaptiveMq(const Rule& r, const ErrorHeuristic& er, std::size_t subdivisions)
        : rule(r), error_heuristic(er), subdivisions(subdivisions) {}

    template<std::size_t DIMBINS, typename F, typename Float, std::size_t DIM, typename Logger>
    auto generate(const std::array<std::size_t, DIMBINS>& bin_resolution,
                  const F& f, const Range<Float, DIM>& range, Logger& logger) const {

        auto r_init = region(f, rule, range.min(), range.max());
        auto errdim_init = error_heuristic(r_init);
        using ERegion = ExtendedRegion<decltype(r_init), decltype(errdim_init)>;

        unsigned int num_threads = std::thread::hardware_concurrency();
        if (num_threads == 0) num_threads = 2;

        Multiqueue<ERegion, MqComparator> mq(4 * num_threads, 2, MqComparator());
        mq.push(ERegion(r_init, errdim_init));

        // Contador atómico compartido entre todos los hilos
        std::atomic<std::size_t> completed{0};  
        std::atomic<bool> stop{false};

        logger.log_progress(std::size_t(0), subdivisions);

        auto worker = [&]() {
            while (true) {
                // Reservar un slot de trabajo
                std::size_t my_slot = completed.fetch_add(1, std::memory_order_relaxed);
                if (my_slot >= subdivisions) {
                    completed.fetch_sub(1, std::memory_order_relaxed);
                    break;
                }

                // Esperar hasta que haya región disponible
                std::optional<ERegion> opt_r;
                while (!opt_r && !stop.load(std::memory_order_relaxed)) {
                    opt_r = mq.try_pop();
                }

                if (!opt_r) break;  // stop fue activado

                auto subregions = opt_r->split(f, std::get<1>(opt_r->extra()));
                for (auto& sr : subregions)
                    mq.push(ERegion(sr, error_heuristic(sr)));

                logger.log_progress(my_slot, subdivisions);
            }
        };

        // Lanzar hilos
        std::vector<std::thread> threads;
        threads.reserve(num_threads);
        for (unsigned int t = 0; t < num_threads; ++t)
            threads.emplace_back(worker);
        for (auto& t : threads)
            t.join();

        // Volcamos la multiqueue al vector de resultado sin ordenar, no modificamos el contenido de mq
        std::vector<ERegion> final_heap;
        final_heap = mq.drain();

        logger.log_progress(subdivisions, subdivisions);
        return final_heap;
    }
};

template<typename Rule, typename ErrorHeuristic>
auto regions_generator_adaptive_multiqueue(const Rule& r, const ErrorHeuristic& er,
                                           std::size_t subdivisions) {
    return RegionsGeneratorAdaptiveMq<Rule, ErrorHeuristic>(r, er, subdivisions);
}

} // namespace viltrum