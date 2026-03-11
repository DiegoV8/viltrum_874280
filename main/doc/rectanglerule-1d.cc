#include "../../viltrum.h"
#include <iostream>
#include <iomanip>

int main() {
    std::array<std::size_t, 1> resolution = {64};

    float sol = viltrum::integrate(
        viltrum::rectangle_rule(resolution),           // Técnica: Regla del Rectángulo
        [](float x) -> float { return std::sin(x); },  // Función: sin(x)
        viltrum::range(0.0f, 3.14159265f)              // Rango: [0, pi]
    );

    std::cout << std::fixed << std::setprecision(6) << sol 
              << " should be much closer to 2 than MC with 64 samples\n";

    return 0;
}