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

        // 1. Calculamos el factor de resolución (igual que en monte-carlo.h)
        double resolution_factor = 1.0;
        for (std::size_t i = 0; i < DIMBINS; ++i) resolution_factor *= bin_resolution[i];

        // 2. El factor final: (Resolución * Volumen del rango) / Muestras totales
        double n_samples = std::pow(steps, DIM);
        double factor = (resolution_factor * range.volume()) / n_samples;

        iterate_grid<DIM>(range, steps, [&](const std::array<Float, DIM>& sample) {
            std::array<std::size_t, DIMBINS> pos;
            for (std::size_t i = 0; i < DIMBINS; ++i) {
                // Si el rango es pequeño (un píxel), la posición es 0 (local al bin)
                // Si el rango es grande (toda la imagen), calculamos el píxel
                if (range.max(i) - range.min(i) <= 1.0001) { 
                    pos[i] = 0; 
                } else {
                    double norm = (sample[i] - range.min(i)) / (range.max(i) - range.min(i));
                    pos[i] = std::min(static_cast<std::size_t>(norm * bin_resolution[i]), bin_resolution[i] - 1);
                }
            }
            
            // IMPORTANTE: El factor de escala en PBRT con IntegratorPerBin
            // debe ser 1.0 / muestras_del_pixel. 
            // Viltrum ya multiplica por 'factor' (resolución total) fuera.
            double n_samples = std::pow(steps, DIM);
            double factor = 1.0 / n_samples;

            bins(pos) += f(sample) * factor;
        });
    }

private:
    template<std::size_t DIM, typename Float, typename Callback>
    void iterate_grid(const Range<Float, DIM>& range, std::size_t steps, Callback cb) const {
        std::array<std::size_t, DIM> indices{};
        std::array<Float, DIM> sample;
        while (true) {
            for (std::size_t d = 0; d < DIM; ++d) {
                Float step_size = (range.max(d) - range.min(d)) / static_cast<Float>(steps);
                sample[d] = range.min(d) + (static_cast<Float>(indices[d]) + 0.5f) * step_size;
            }
            cb(sample);
            std::size_t i = 0;
            while (true) {
                if (++indices[i] < steps) break;
                indices[i] = 0;
                if (++i == DIM) return;
            }
        }
    }
};

inline auto rectangle_rule(std::size_t steps) { return RectangleRule(steps); }

}