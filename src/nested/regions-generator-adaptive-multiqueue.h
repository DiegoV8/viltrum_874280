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
        auto errA = std::get<0>(a.extra());
        auto errB = std::get<0>(b.extra());
        
        bool nanA = std::isnan(errA);
        bool nanB = std::isnan(errB);
        
        // Si ambos son NaN, son equivalentes (ninguno es menor que el otro)
        if (nanA && nanB) return false; 
        
        // Si solo uno es NaN, lo mandamos al fondo (menor prioridad)
        if (nanA) return true;  
        if (nanB) return false;
        
        return errA < errB;
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

        #pragma omp parallel
        {
            while (true) {
                // ¿Queda trabajo por hacer?
                std::size_t current = completed.load(std::memory_order_relaxed);
                if (current >= subdivisions) break;

                // Intentar obtener una región de forma segura de nuestra Multiqueue
                std::optional<ERegion> opt_r = mq.try_pop();
                if (!opt_r) {
                    // No hay nada ahora, cedemos el turno brevemente
                    std::this_thread::yield();
                    continue;
                }

                // Reservar el slot DESPUÉS de tener trabajo real
                std::size_t my_slot = completed.fetch_add(1, std::memory_order_relaxed);
                if (my_slot >= subdivisions) {
                    // Nos pasamos, devolver la región y salir
                    mq.push(*opt_r);
                    break;
                }

                // ¡Aquí f se ejecuta en paralelo de forma 100% segura bajo OpenMP!
                auto subregions = opt_r->split(f, std::get<1>(opt_r->extra()));
                for (auto& sr : subregions)
                    mq.push(ERegion(sr, error_heuristic(sr)));

                #pragma omp critical(logger_update)
                {
                    logger.log_progress(my_slot + 1, subdivisions);
                }
            }
        } // Fin de la región paralela de OpenMP (barrera implícita, todos los hilos coordinados aquí)

        // Volcamos la multiqueue al vector de resultado sin ordenar
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