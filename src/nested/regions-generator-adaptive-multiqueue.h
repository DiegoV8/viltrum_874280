#pragma once
#include "nested.h"
#include "error-heuristic.h"
#include <thread>
#include "../relaxedPriorityQueues_874280/src/multiqueue/multiqueue.hpp"

struct MqComparator {
    template<typename T>
    bool operator()(const T& a, const T& b) const {
        return std::get<0>(a.extra()) < std::get<0>(b.extra());
    }
};

namespace viltrum {

template<typename Rule, typename ErrorHeuristic, typename = std::enable_if_t<is_nested<Rule>::value>>
class RegionsGeneratorAdaptiveMq {
    Rule rule;
    ErrorHeuristic error_heuristic;
    std::size_t subdivisions;

public:
    RegionsGeneratorAdaptiveMq(const Rule& r, const ErrorHeuristic& er, std::size_t subdivisions) 
        : rule(r), error_heuristic(er), subdivisions(subdivisions) {}

    template<std::size_t DIMBINS, typename F, typename Float, std::size_t DIM, typename Logger>
    auto generate(const std::array<std::size_t,DIMBINS>& bin_resolution,
		const F& f, const Range<Float,DIM>& range, Logger& logger) const {

        auto r_init = region(f, rule, range.min(), range.max());
        auto errdim_init = error_heuristic(r_init);
        using ERegion = ExtendedRegion<decltype(r_init), decltype(errdim_init)>;

        // Numero de colas depende del numero de threads, c=2
        unsigned int num_threads = std::thread::hardware_concurrency();
        if (num_threads == 0) num_threads = 2;
        Multiqueue<ERegion, MqComparator> mq(2*num_threads, 2, MqComparator());

        logger.log_progress(std::size_t(0), subdivisions);
        mq.push(ERegion(r_init, errdim_init));

        for (std::size_t i = 0; i < subdivisions; ++i) {
            auto opt_r = mq.try_pop();
            if (!opt_r) break;

            ERegion current_r = *opt_r;
            auto subregions = current_r.split(f, std::get<1>(current_r.extra()));
            for (auto& sr : subregions) {
                auto errdim = error_heuristic(sr);
                mq.push(ERegion(sr, errdim));
            }
            
            logger.log_progress(i, subdivisions);
        }

        // Reconvertir la Multiqueue a un vector
        std::vector<ERegion> final_heap;
        while (auto final_r = mq.try_pop()) {
            final_heap.push_back(*final_r);
        }

        logger.log_progress(subdivisions, subdivisions);
        return final_heap;
	    }
};

template<typename Rule, typename ErrorHeuristic>
auto regions_generator_adaptive_multiqueue(const Rule& r, const ErrorHeuristic& er, std::size_t subdivisions) {
    return RegionsGeneratorAdaptiveMq<Rule,ErrorHeuristic>(r,er,subdivisions);
}


}