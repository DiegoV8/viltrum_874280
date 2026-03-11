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

    // Versión para salida única
    template<typename F, typename Float, typename Logger>
    auto integrate(const F& f, const Range<Float, DIM>& range, Logger& logger) const {
        using ReturnType = decltype(f(std::declval<std::array<Float, DIM>>()));
        ReturnType total_sum = ReturnType(0);
        
        std::array<std::size_t, 0> no_res;
        auto single_bin = [&](const std::array<std::size_t, 0>&) -> ReturnType& { return total_sum; };
        
        integrate_impl<decltype(single_bin), 0, F, Float, Logger>(single_bin, no_res, f, range, logger);
        return total_sum;
    }

    // Versión para Bins (Corregida para que coincida con lo que busca viltrum::integrate)
    template<typename Bins, std::size_t DIMBINS, typename F, typename Float, typename Logger>
    void integrate(Bins& bins, const std::array<std::size_t, DIMBINS>& bin_res, const F& f, const Range<Float, DIM>& range, Logger& logger) const {
        integrate_impl<Bins, DIMBINS, F, Float, Logger>(bins, bin_res, f, range, logger);
    }

private:
    template<typename Bins, std::size_t DIMBINS, typename F, typename Float, typename Logger>
    void integrate_impl(Bins& bins, const std::array<std::size_t, DIMBINS>& bin_resolution, const F& f, const Range<Float, DIM>& range, Logger& logger) const {
        using FloatType = Float;
        
        std::array<FloatType, DIM> step;
        FloatType cell_volume = 1.0;
        unsigned long total_cells = 1;

        for (std::size_t d = 0; d < DIM; ++d) {
            step[d] = (range.max(d) - range.min(d)) / FloatType(resolution[d]);
            cell_volume *= step[d];
            total_cells *= resolution[d];
        }

        for (unsigned long i = 0; i < total_cells; ++i) {
            logger.log_progress(i, total_cells);
            std::array<FloatType, DIM> sample_point;
            unsigned long temp_i = i;

            for (std::size_t d = 0; d < DIM; ++d) {
                std::size_t grid_pos = temp_i % resolution[d];
                temp_i /= resolution[d];
                sample_point[d] = range.min(d) + (FloatType(grid_pos) + 0.5f) * step[d];
            }

            if constexpr (DIMBINS > 0) {
                std::array<std::size_t, DIMBINS> bin_pos;
                for (std::size_t d = 0; d < DIMBINS; ++d) {
                    // Mapeo del punto de la regla al bin correspondiente
                    bin_pos[d] = std::min(std::size_t(bin_resolution[d] * (sample_point[d] - range.min(d)) / (range.max(d) - range.min(d))), bin_resolution[d] - 1);
                }
                bins(bin_pos) += f(sample_point) * cell_volume;
            } else {
                std::array<std::size_t, 0> dummy;
                bins(dummy) += f(sample_point) * cell_volume;
            }
        }
    }
};

// Helper para crear la regla
template<std::size_t DIM>
RectangleRule<DIM> rectangle_rule(std::array<std::size_t, DIM> res) {
    return RectangleRule<DIM>(res);
}

} // namespace viltrum