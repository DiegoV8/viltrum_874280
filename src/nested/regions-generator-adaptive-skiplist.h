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
        // Para que la Skiplist se comporte como una cola de prioridad de máximos
        // (extrayendo el mayor error primero desde la cabeza), usamos '>'
        // Si se usa '<', se extraería el error más pequeño.
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

        // 1. Inicialización de la región base
        auto r_init = region(f, rule, range.min(), range.max());
        auto errdim_init = error_heuristic(r_init);
        
        // Definimos el tipo de región extendida y el puntero inteligente para la Skiplist
        using ERegion = ExtendedRegion<decltype(r_init), decltype(errdim_init)>;
        using ERegionPtr = std::shared_ptr<ERegion>;

        // 2. Configuración de concurrencia
        unsigned int num_threads = std::thread::hardware_concurrency();
        if (num_threads == 0) num_threads = 2;

        // 3. Inicialización de la Skiplist (usando {} para evitar el "most vexing parse")
        // Se recomienda un factor de relajación para mejorar el rendimiento en entornos concurrentes
        Skiplist<ERegionPtr, SlComparator> sl{SlComparator(), int(num_threads)*128};
        sl.push(std::make_shared<ERegion>(r_init, errdim_init));

        // Contador atómico para controlar el número de subdivisiones realizadas
        std::atomic<std::size_t> done{0};
        logger.log_progress(std::size_t(0), subdivisions);

        // 4. Definición del Worker para procesamiento paralelo
        auto worker = [&]() {
            while (true) {
                // Reservamos un slot de trabajo
                std::size_t my_slot = done.fetch_add(1, std::memory_order_relaxed);
                if (my_slot >= subdivisions) {
                    done.fetch_sub(1, std::memory_order_relaxed);
                    break;
                }

                // Extracción de la región con mayor error (espera activa si es necesario)
                std::optional<ERegionPtr> opt_ptr;
                while (!opt_ptr) {
                    opt_ptr = sl.try_pop(); //
                    // Si la cola está vacía pero el trabajo global no ha terminado, reintentamos.
                    // Si ya se alcanzó el límite de subdivisiones, salimos.
                    if (!opt_ptr && done.load() >= subdivisions) break;
                }

                if (!opt_ptr) {
                    done.fetch_sub(1, std::memory_order_relaxed);
                    break;
                }

                // 5. División de la región y actualización de la estructura
                auto subregions = (*opt_ptr)->split(f, std::get<1>((*opt_ptr)->extra()));
                for (auto& sr : subregions) {
                    sl.push(std::make_shared<ERegion>(sr, error_heuristic(sr))); //
                }

                logger.log_progress(my_slot, subdivisions);
            }
        };

        // 6. Lanzamiento y sincronización de hilos
        std::vector<std::thread> threads;
        threads.reserve(num_threads);
        for (unsigned int t = 0; t < num_threads; ++t)
            threads.emplace_back(worker);
        
        for (auto& t : threads)
            t.join();

        // 7. Recolección final de resultados
        // Vaciamos la skiplist y convertimos los shared_ptr de vuelta a objetos ERegion
        auto final_ptrs = sl.drain(); //
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

} // namespace viltrum