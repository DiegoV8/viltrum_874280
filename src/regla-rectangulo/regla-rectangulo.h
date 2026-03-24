#pragma once
#include <array>
#include <cmath>
#include <vector>
#include <algorithm>
#include <limits>
#include <random>
#include "../range.h"

namespace viltrum {

template<typename Float, std::size_t DIM>
struct RectangleSequence {
    std::array<Float, DIM> data;

    // Magia de PBRT: [0] y [1] son la cámara. [2+] son los rebotes de luz.
    Float operator[](std::size_t i) const {
    // Si la dimensión está dentro de las calculadas por la rejilla (ej. x, y del píxel)
    if (i < DIM) {
        return data[i]; 
    }
    
    // Para dimensiones extra (rebotes, luces, etc.), usamos el centro (0.5)
    // Esto equivale a evaluar el "hiper-centro" de la caja en esas dimensiones.
    return static_cast<Float>(0.5);
}

    struct Iterator {
        const RectangleSequence* seq;
        std::size_t index;

        Iterator(const RectangleSequence* s, std::size_t i) : seq(s), index(i) {}

        Float operator*() const { return (*seq)[index]; }
        Iterator& operator++() { ++index; return *this; }
        Iterator operator++(int) { Iterator tmp = *this; ++(*this); return tmp; }
        bool operator==(const Iterator& o) const { return index == o.index; }
        bool operator!=(const Iterator& o) const { return index != o.index; }
    };

    Iterator begin() const { return Iterator{this, 0}; }
    Iterator end() const { return Iterator{this, std::numeric_limits<std::size_t>::max()}; }
};

class RectangleRule {
    std::size_t steps;

public:
    RectangleRule(std::size_t s) : steps(s) {}

    template<typename Bins, std::size_t DIMBINS, typename F, typename Float, std::size_t DIM, typename Logger>
    void integrate(Bins& bins, const std::array<std::size_t, DIMBINS>& bin_resolution,
                   const F& f, const Range<Float, DIM>& range, Logger& logger) const {
        
        if (steps == 0) return;

        // 1. Número total de muestras en la rejilla
        double n_samples = std::pow(steps, DIM);

        // 2. Calculamos el factor de resolución
        double resolution_factor = 1.0;
        for (std::size_t i = 0; i < DIMBINS; ++i) {
            resolution_factor *= bin_resolution[i];
        }

        // 3. Factor final (escala de PBRT y Viltrum)
        double factor = (resolution_factor * range.volume()) / n_samples;
        
        unsigned long current_sample = 0;

        // 4. Llamamos a la función recursiva
        iterate_recursive<Float, DIM>(range, steps, [&](const std::array<Float, DIM>& sample) {
            logger.log_progress(current_sample++, static_cast<unsigned long>(n_samples));

            if (range.is_inside(sample)) {
                std::array<std::size_t, DIMBINS> pos;
                for (std::size_t i = 0; i < DIMBINS; ++i) {
                    pos[i] = std::size_t(bin_resolution[i] * (sample[i] - range.min(i)) / (range.max(i) - range.min(i)));
                }
                
                // ---- MODIFICACIÓN AQUÍ ----
                // Hemos eliminado RectangleSequence. 
                // Pasamos 'sample' (el std::array crudo) directamente a 'f'.
                // Fubini pasará este array al Monte Carlo anidado, y será Monte Carlo 
                // quien genere la secuencia infinita final para PBRT.
                bins(pos) += f(sample) * factor;
            }
        });

        logger.log_progress(static_cast<unsigned long>(n_samples), static_cast<unsigned long>(n_samples));
    }

private:
    // La función auxiliar que hace la magia de la rejilla multidimensional
    template<typename Float, std::size_t DIM, typename Callback>
    void iterate_recursive(const Range<Float, DIM>& range, std::size_t steps_per_dim, 
                           Callback cb, std::array<Float, DIM> current_pos = {}, std::size_t depth = 0) const {
        if (depth == DIM) {
            cb(current_pos);
            return;
        }

        Float step_size = (range.max(depth) - range.min(depth)) / steps_per_dim;
        for (std::size_t i = 0; i < steps_per_dim; ++i) {
            // Evaluamos en el punto medio de cada rectángulo
            current_pos[depth] = range.min(depth) + (i + 0.5f) * step_size;
            iterate_recursive<Float, DIM>(range, steps_per_dim, cb, current_pos, depth + 1);
        }
    }
};

inline auto rectangle_rule(std::size_t steps) { return RectangleRule(steps); }

}