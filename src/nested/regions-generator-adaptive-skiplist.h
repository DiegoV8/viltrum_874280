#pragma once
#include "nested.h"
#include "error-heuristic.h"
#include <memory>
#include <thread>
#include <utility>
#include <atomic>
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

        // 1. Inicialización de la primera región
        auto r_init = region(f, rule, range.min(), range.max());
        auto errdim_init = error_heuristic(r_init);
        using ERegion  = ExtendedRegion<decltype(r_init), decltype(errdim_init)>;
        using PERegion = std::shared_ptr<ERegion>;

        unsigned int num_threads = std::thread::hardware_concurrency();
        if (num_threads == 0) num_threads = 2;

        // Skiplist con factor de relajación para evitar contención
        Skiplist<PERegion, SlComparator> sl(SlComparator(), num_threads * 8);
        sl.push(std::make_shared<ERegion>(r_init, errdim_init));

        // 2. Control de terminación concurrente
        std::atomic<std::size_t> completed_subdivisions{0};
        std::atomic<int> active_workers{0};

        auto worker_task = [&]() {
            while (true) {
                // Condición de salida global: ¿se ha alcanzado el objetivo?
                if (completed_subdivisions.load(std::memory_order_relaxed) >= subdivisions) break;

                auto opt_r = sl.try_pop();

                if (!opt_r) {
                    // Si la lista está vacía, verificamos si otros hilos aún están procesando.
                    // Si nadie está trabajando (active_workers == 0), no habrá nuevos datos.
                    if (active_workers.load(std::memory_order_acquire) == 0 && sl.empty()) {
                        break; 
                    }
                    std::this_thread::yield(); // Espera activa ligera
                    continue;
                }

                // Marcamos que este hilo está produciendo nuevas regiones
                active_workers.fetch_add(1, std::memory_order_acq_rel);

                // Doble comprobación atómica para no exceder las subdivisiones
                if (completed_subdivisions.fetch_add(1, std::memory_order_acq_rel) >= subdivisions) {
                    sl.push(*opt_r); // Devolvemos la región si nos pasamos
                    completed_subdivisions.fetch_sub(1, std::memory_order_relaxed);
                    active_workers.fetch_sub(1, std::memory_order_release);
                    break;
                }

                // 3. Procesamiento de la región (Split)
                ERegion& current_r = **opt_r;
                auto subregions = current_r.split(f, std::get<1>(current_r.extra()));

                for (auto& sr : subregions) {
                    auto errdim = error_heuristic(sr);
                    sl.push(std::make_shared<ERegion>(std::move(sr), std::move(errdim)));
                }

                active_workers.fetch_sub(1, std::memory_order_release);
                
                // Log de progreso solo en un hilo o con muestreo para no saturar la salida
                if (completed_subdivisions.load() % 10 == 0) {
                    logger.log_progress(completed_subdivisions.load(), subdivisions);
                }
            }
        };

        // 4. Lanzamiento de hilos
        std::vector<std::thread> threads;
        for (unsigned int t = 0; t < num_threads; ++t) {
            threads.emplace_back(worker_task);
        }

        for (auto& t : threads) t.join();

        // 5. Recolección final mediante drain()
        std::vector<ERegion> final_heap;
        auto drained_ptrs = sl.drain();
        final_heap.reserve(drained_ptrs.size());
        for(auto& ptr : drained_ptrs) {
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

} // namespace viltrum