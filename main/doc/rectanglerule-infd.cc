#include "../../viltrum.h"
#include <iostream>
#include <iomanip>
#include <array>
#include <algorithm>
#include <type_traits>

int main() {
    float decaying_factor = 0.75f;
    
    // Dimensionalidad
    const std::size_t D = 6;

    auto integrand = [decaying_factor] (const std::array<float, D>& seq) -> float {
            auto x = seq.begin();
            float sum = 0.0f;
            float term = 1.0f;
        
        while (x != seq.end() && (*x) < decaying_factor) { 
            float val_rr = *x;
            ++x; 
            if (x == seq.end()) break;
            
            float val_term = *x;
            term *= 2.0f * val_term; 
            sum += term; 
            ++x;
        }
        return sum;
    };

    std::array<std::size_t, D> resolution;
    resolution.fill(4); // 4 divisiones por dimensión = 4^6 muestras totales

    // efinimos los límites del hipercubo unitario [0, 1]^D
    std::array<float, D> min_bounds;
    std::array<float, D> max_bounds;
    min_bounds.fill(0.0f);
    max_bounds.fill(1.0f);

    float sol = viltrum::integrate(
        viltrum::rectangle_rule(resolution), 
        integrand,
        viltrum::range(min_bounds, max_bounds) 
    );

    std::cout << std::fixed << std::setprecision(6);
    std::cout << "Integral (Rect Rule D=" << D << "): " << sol << "\n";
    std::cout << "Target Value: " << (decaying_factor / (1.0f - decaying_factor)) << std::endl;

    return 0;
}

