#include "../../viltrum.h"
#include <iostream>
#include <vector>

int main() {
    // Definimos el integrando f(x,y) = x^2 + y^2
    auto integrand = [] (const std::array<float,2>& x) -> float { 
        return x[0]*x[0] + x[1]*x[1]; 
    };

    auto range = viltrum::range(std::array<float,2>{0.0f, 0.0f}, std::array<float,2>{1.0f, 1.0f}); 

    std::vector<float> sol_vec(10, 0.0f);

    viltrum::integrate(
        viltrum::rectangle_rule(64), 
        sol_vec, 
        integrand,
        range
    );

    for (std::size_t i = 0; i < sol_vec.size(); ++i) { 
        std::cout << "Bin " << i << ": " << sol_vec[i] << "\n";
    }

    return 0;
}