#include "../../viltrum.h"
#include <iostream>
#include <vector>
#include <array>

int main() {
    // 1. Definimos los contenedores para los resultados (bins)
    float sol[10] = {0.0f}; 
    std::vector<float> sol_vec(10, 0.0f);

    // 2. Definimos el acceso a los bins (mapeo de posición a memoria)
    auto sol_access = [&sol] (const std::array<std::size_t, 1>& pos) -> float& { 
        return sol[pos[0]]; 
    }; 

    // 3. Función e Integrando (2D) y Rango
    auto integrand = [] (const std::array<float, 2>& x) -> float { 
        return x[0]*x[0] + x[1]*x[1]; 
    };
    auto range = viltrum::range(std::array<float, 2>{0.0f, 0.0f}, std::array<float, 2>{1.0f, 1.0f});

    // 4. Configuración de la Regla del Rectángulo
    // Definimos una rejilla de, por ejemplo, 100x100 para tener 10.000 muestras deterministas
    std::array<std::size_t, 2> grid_resolution = {100, 100};
    auto rect_rule = viltrum::rectangle_rule(grid_resolution);

    // 5. Integración con acceso explícito a bins
    viltrum::integrate(
        rect_rule,
        sol_access,
        std::array<std::size_t, 1>{10}, // Resolución de los bins (1D)
        integrand,
        range
    );

    // Salida de resultados
    for (std::size_t i = 0; i < 10; ++i) { 
        std::cout << "Bin " << i << ": " << sol[i] << "\n";
    }

    // 6. Integración con contenedor implícito (std::vector)
    // viltrum::integrate debe estar sobrecargado para manejar vectores
    viltrum::integrate(
        rect_rule,
        sol_vec, 
        integrand,
        range
    );

    return 0;
}