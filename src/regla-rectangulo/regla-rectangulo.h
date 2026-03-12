#pragma once
#include <array>
#include <cmath>
#include <vector>
#include "../range.h"

namespace viltrum {

class RectangleRule {
    std::size_t steps;

public:
    RectangleRule(std::size_t s) : steps(s) {}

    template<typename Bins, std::size_t DIMBINS, typename F, typename Float, std::size_t DIM, typename Logger>
    void integrate(Bins& bins, const std::array<std::size_t, DIMBINS>& bin_resolution,
        const F& f, const Range<Float, DIM>& range, Logger& logger) const {

        // 1. Evitar división por cero si steps es 0 (por seguridad)
        if (steps == 0) return;

        double step_volume = range.volume() / std::pow(static_cast<double>(steps), static_cast<double>(DIM));
        
        double resolution_factor = 1;
        for (std::size_t i = 0; i < DIMBINS; ++i) resolution_factor *= bin_resolution[i];
        double factor = resolution_factor * step_volume;

        // Log inicial
        logger.log_progress(0, 1);

        iterate_grid<DIM>(range, steps, [&](const std::array<Float, DIM>& sample) {
            // El integrando de PBRT (f) espera una secuencia. 
            // Pasamos 'sample' que es un std::array (cumple con ser una clase/contenedor).
            if (range.is_inside(sample)) {
                std::array<std::size_t, DIMBINS> pos;
                for (std::size_t i = 0; i < DIMBINS; ++i) {
                    // 2. Clamp preventivo para evitar desbordamiento de índice por precisión flotante
                    double norm = (sample[i] - range.min(i)) / (range.max(i) - range.min(i));
                    pos[i] = std::min(static_cast<std::size_t>(bin_resolution[i] * norm), bin_resolution[i] - 1);
                }
                
                // f(sample) llama al trazador de rayos. 
                // IMPORTANTE: Asegúrate de que DIM en el main coincida con lo que espera f.
                bins(pos) += f(sample) * factor;
            }
        });

        logger.log_progress(1, 1);
    }

private:
    template<std::size_t DIM, typename Float, typename Callback>
    void iterate_grid(const Range<Float, DIM>& range, std::size_t steps, Callback cb) const {
        if (steps == 0) return;
        std::array<std::size_t, DIM> indices{};
        std::array<Float, DIM> sample;
        
        while (true) {
            for (std::size_t d = 0; d < DIM; ++d) {
                Float step_size = (range.max(d) - range.min(d)) / static_cast<Float>(steps);
                // Punto medio para mayor precisión
                sample[d] = range.min(d) + (static_cast<Float>(indices[d]) + 0.5f) * step_size;
            }
            cb(sample);

            std::size_t i = 0;
            // Lógica de acarreo para simular bucles anidados
            while (++indices[i] == steps) {
                indices[i] = 0;
                if (++i == DIM) return; 
            }
        }
    }
};

inline auto rectangle_rule(std::size_t steps) {
    return RectangleRule(steps);
}

}