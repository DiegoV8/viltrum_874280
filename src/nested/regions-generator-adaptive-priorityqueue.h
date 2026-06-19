#pragma once
#include "nested.h"
#include "error-heuristic.h"
#include <thread>
#include <atomic>
#include <vector>
#include "../relaxedPriorityQueues_874280/src/multiqueue/multiqueue.hpp"

// Estructura usada para comparar dos elementos en la estructura.
struct PqComparator {
    template<typename T>
    bool operator()(const T& a, const T& b) const {
        auto errA = std::get<0>(a.extra());
        auto errB = std::get<0>(b.extra());
        
        bool nanA = std::isnan(errA);
        bool nanB = std::isnan(errB);
        
        if (nanA && nanB) return false; 
        
        if (nanA) return true;  
        if (nanB) return false;
        
        return errA < errB;
    }
};

namespace viltrum {

template<typename Rule, typename ErrorHeuristic,
         typename = std::enable_if_t<is_nested<Rule>::value>>
class RegionsGeneratorAdaptivePq {
    Rule rule;
    ErrorHeuristic error_heuristic;
    std::size_t subdivisions;

public:
    RegionsGeneratorAdaptivePq(const Rule& r, const ErrorHeuristic& er, std::size_t subdivisions)
        : rule(r), error_heuristic(er), subdivisions(subdivisions) {}

    template<std::size_t DIMBINS, typename F, typename Float, std::size_t DIM, typename Logger>
    auto generate(const std::array<std::size_t, DIMBINS>& bin_resolution,
                  const F& f, const Range<Float, DIM>& range, Logger& logger) const {

        auto r_init = region(f, rule, range.min(), range.max());
        auto errdim_init = error_heuristic(r_init);
        using ERegion = ExtendedRegion<decltype(r_init), decltype(errdim_init)>;

        unsigned int num_threads = std::thread::hardware_concurrency();
        if (num_threads == 0) num_threads = 2;

        Multiqueue<ERegion, PqComparator> mq(1, 1, PqComparator()); // Equivalente a una cola
        mq.push(ERegion(r_init, errdim_init));

        std::atomic<std::size_t> completed{0};  
        std::atomic<bool> stop{false};

        logger.log_progress(std::size_t(0), subdivisions);

        #pragma omp parallel
        {
            while (true) {
                std::size_t current = completed.load(std::memory_order_relaxed);
                if (current >= subdivisions) break;

                std::optional<ERegion> opt_r = mq.try_pop();
                if (!opt_r) {
                    std::this_thread::yield();
                    continue;
                }

                std::size_t my_slot = completed.fetch_add(1, std::memory_order_relaxed);
                if (my_slot >= subdivisions) {
                    mq.push(*opt_r);
                    break;
                }

                auto subregions = opt_r->split(f, std::get<1>(opt_r->extra()));
                for (auto& sr : subregions)
                    mq.push(ERegion(sr, error_heuristic(sr)));

                #pragma omp critical(logger_update)
                {
                    logger.log_progress(my_slot + 1, subdivisions);
                }
            }
        }

        std::vector<ERegion> final_heap;
        final_heap = mq.drain();

        logger.log_progress(subdivisions, subdivisions);

        return final_heap;
    }
};

template<typename Rule, typename ErrorHeuristic>
auto regions_generator_adaptive_priorityqueue(const Rule& r, const ErrorHeuristic& er,
                                           std::size_t subdivisions) {
    return RegionsGeneratorAdaptiveMq<Rule, ErrorHeuristic>(r, er, subdivisions);
}

}