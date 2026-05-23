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

        // Contador atómico compartido entre todos los hilos
        std::atomic<std::size_t> completed{0};

        logger.log_progress(std::size_t(0), subdivisions);

        // 4. Región paralela OpenMP — análoga al multiqueue
        #pragma omp parallel
        {
            while (true) {
                // ¿Queda trabajo por hacer?
                std::size_t current = completed.load(std::memory_order_relaxed);
                if (current >= subdivisions) break;

                // Intentar obtener una región de forma segura de la Skiplist
                std::optional<ERegionPtr> opt_ptr = sl.try_pop();
                if (!opt_ptr) {
                    // No hay nada ahora, cedemos el turno brevemente
                    std::this_thread::yield();
                    continue;
                }

                // Reservar el slot DESPUÉS de tener trabajo real (evita la race condition
                // del diseño original que reservaba antes de hacer el pop)
                std::size_t my_slot = completed.fetch_add(1, std::memory_order_relaxed);
                if (my_slot >= subdivisions) {
                    // Nos pasamos: devolver la región y salir
                    sl.push(*opt_ptr);
                    break;
                }

                // División de la región y actualización de la estructura
                // f se ejecuta en paralelo de forma segura bajo OpenMP
                auto subregions = (*opt_ptr)->split(f, std::get<1>((*opt_ptr)->extra()));
                for (auto& sr : subregions) {
                    sl.push(std::make_shared<ERegion>(sr, error_heuristic(sr)));
                }

                #pragma omp critical(logger_update)
                {
                    logger.log_progress(my_slot + 1, subdivisions);
                }
            }
        } // Fin de la región paralela de OpenMP (barrera implícita)

        // 5. Recolección final de resultados
        // Vaciamos la skiplist y convertimos los shared_ptr de vuelta a objetos ERegion
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

} // namespace viltrum