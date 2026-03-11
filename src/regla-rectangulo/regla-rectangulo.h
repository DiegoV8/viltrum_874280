#pragma once
#include <array>
#include <vector>
#include <algorithm>
#include "../range.h"

namespace viltrum {

template<std::size_t DIM>
class RectangleRule {
    std::array<std::size_t, DIM> resolution;

public:
    RectangleRule(const std::array<std::size_t, DIM>& res) : resolution(res) {}

    // Esta es la función que viltrum::integrate buscará
    template<typename Bins, std::size_t DIMBINS, typename F, typename Float, typename Logger>
    void integrate(Bins& bins, const std::array<std::size_t, DIMBINS>& bin_resolution, 
                   const F& f, const Range<Float, DIM>& range, Logger& logger) const {
        
        std::array<Float, DIM> step;
        Float cell_volume = 1.0;
        unsigned long total_cells = 1;

        for (std::size_t i = 0; i < DIM; ++i) {
            step[i] = (range.max(i) - range.min(i)) / Float(resolution[i]);
            cell_volume *= step[i];
            total_cells *= resolution[i];
        }

        // El factor de normalización para que la energía sea correcta
        // En PBRT, si lanzas N muestras, cada una debe valer 1/N
        Float total_factor = cell_volume; 

        for (unsigned long i = 0; i < total_cells; ++i) {
            logger.log_progress(i, total_cells);
            
            std::array<Float, DIM> sample_point;
            unsigned long temp_i = i;

            for (std::size_t d = 0; d < DIM; ++d) {
                std::size_t grid_pos = temp_i % resolution[d];
                temp_i /= resolution[d];
                // Punto medio de la celda de la rejilla
                sample_point[d] = range.min(d) + (Float(grid_pos) + 0.5f) * step[d];
            }

            // Calculamos a qué píxel (bin) corresponde esta muestra
            if constexpr (DIMBINS > 0) {
                std::array<std::size_t, DIMBINS> bin_pos;
                for (std::size_t d = 0; d < DIMBINS; ++d) {
                    Float normalized = (sample_point[d] - range.min(d)) / (range.max(d) - range.min(d));
                    bin_pos[d] = std::min(std::size_t(normalized * bin_resolution[d]), bin_resolution[d] - 1);
                }
                // Llamamos a la lambda bins_accessor del .cc
                bins(bin_pos) += f(sample_point) * total_factor;
            }
        }
    }
};

template<std::size_t DIM>
RectangleRule<DIM> rectangle_rule(std::array<std::size_t, DIM> res) {
    return RectangleRule<DIM>(res);
}

} // namespace viltrum