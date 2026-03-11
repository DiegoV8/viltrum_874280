#pragma once
#include <array>
#include <vector>
#include "../range.h"

namespace viltrum {

template<std::size_t DIM>
class RectangleRule {
    std::array<std::size_t, DIM> resolution;

public:
    RectangleRule(const std::array<std::size_t, DIM>& res) : resolution(res) {}

    /**
     * Versión para salida única (Single Output)
     * Requerida por: auto sol = viltrum::integrate(rect, f, range);
     */
    template<typename F, typename Float, typename Logger>
    Float integrate(const F& f, const Range<Float, DIM>& range, Logger& logger) const {
        Float total_sum = 0;
        
        // Creamos un bin ficticio de tamaño 1 para reutilizar la lógica
        auto single_bin = [&](const std::array<std::size_t, 0>&) -> Float& { return total_sum; };
        std::array<std::size_t, 0> no_res;
        
        // Llamamos a la lógica principal (DIMBINS = 0 indica salida única)
        integrate_impl(single_bin, no_res, f, range, logger);
        
        return total_sum;
    }

    /**
     * Versión para Bins (la que ya tenías, pero como implementación privada/interna)
     */
    template<typename Bins, std::size_t DIMBINS, typename F, typename Float, typename Logger>
    void integrate(Bins& bins, const std::array<std::size_t, DIMBINS>& bin_resolution,
                  const F& f, const Range<Float, DIM>& range, Logger& logger) const {
        integrate_impl(bins, bin_resolution, f, range, logger);
    }

private:
    // Lógica centralizada para evitar duplicar el cálculo de la rejilla
    template<typename Bins, std::size_t DIMBINS, typename F, typename Float, typename Logger>
    void integrate_impl(Bins& bins, const std::array<std::size_t, DIMBINS>& bin_resolution,
                       const F& f, const Range<Float, DIM>& range, Logger& logger) const {
        
        double cell_volume = range.volume();
        unsigned long total_cells = 1;
        std::array<Float, DIM> step;

        for (std::size_t i = 0; i < DIM; ++i) {
            step[i] = (range.max(i) - range.min(i)) / Float(resolution[i]);
            cell_volume /= double(resolution[i]);
            total_cells *= resolution[i];
        }

        double bin_factor = 1.0;
        for (std::size_t i = 0; i < DIMBINS; ++i) bin_factor *= bin_resolution[i];
        double total_factor = bin_factor * cell_volume;

        for (unsigned long i = 0; i < total_cells; ++i) {
            logger.log_progress(i, total_cells);
            std::array<Float, DIM> sample_point;
            unsigned long temp_i = i;

            for (std::size_t d = 0; d < DIM; ++d) {
                std::size_t grid_pos = temp_i % resolution[d];
                temp_i /= resolution[d];
                sample_point[d] = range.min(d) + (Float(grid_pos) + 0.5f) * step[d];
            }

            // Invocación: f siempre recibe el std::array (la secuencia)
            if (DIMBINS > 0) {
                std::array<std::size_t, DIMBINS> bin_pos;
                for (std::size_t d = 0; d < DIMBINS; ++d) {
                    bin_pos[d] = std::size_t(bin_resolution[d] * (sample_point[d] - range.min(d)) / (range.max(d) - range.min(d)));
                }
                bins(bin_pos) += f(sample_point) * total_factor;
            } else {
                // Para single output, el factor de bin es 1
                bins({}) += f(sample_point) * cell_volume;
            }
        }
        logger.log_progress(total_cells, total_cells);
    }
};

// Helper
template<std::size_t DIM>
auto rectangle_rule(const std::array<std::size_t, DIM>& resolution) {
    return RectangleRule<DIM>(resolution);
}

}