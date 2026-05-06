#include <iostream>
#include <vector>
#include <chrono>
#include <iomanip>
#include <fstream>
#include <filesystem>

#include "../../viltrum.h"

using namespace viltrum;
using namespace viltrumtest;

// Función genérica para testear cualquier integrando
template<typename Function, std::size_t N>
void run_test(const std::string& name, const Function& f, const std::vector<std::size_t>& steps, std::ofstream& file, const std::string& run_name) {
    auto range = range_all<N>(0.0f, 1.0f);
    float real_value = f.integral_primary();

    std::cout << "\n--- Test: " << name << " ---" << std::endl;
    std::cout << "Valor analítico real: " << real_value << std::endl;
    std::cout << std::left << std::setw(12) << "Iteraciones" 
              << "| " << std::setw(15) << "Resultado" 
              << "| " << std::setw(15) << "Error Rel." 
              << "| " << std::setw(15) << "Tiempo (ms)" << std::endl;
    std::cout << "-----------------------------------------------------------------------" << std::endl;

    for (auto s : steps) {
        auto integrator = integrator_adaptive_iterations_parallel(
            nested(simpson, trapezoidal), 
            s
        );

        // Medición de tiempo de ejecución
        auto start = std::chrono::high_resolution_clock::now();
        float result = integrate(integrator, f, range);
        auto end = std::chrono::high_resolution_clock::now();
        
        std::chrono::duration<double, std::milli> duration = end - start;

        float error = std::abs(result - real_value) / (real_value != 0 ? real_value : 1.0f);

        if (file.is_open()) {
            file << run_name << "_" << name << ";" 
                    << s << ";"               
                    << error << ";" 
                    << duration.count() << "\n";
        }

        std::cout << std::left << std::setw(12) << s 
                  << "| " << std::setw(15) << result 
                  << "| " << std::setw(15) << error 
                  << "| " << std::setw(15) << duration.count() << std::endl;
    }
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Uso: " << argv[0] << " <nombre_del_test>" << std::endl;
        return 1;
    }
    std::string run_name = argv[1];
    std::string path = "../main/function-tests/resultados/res.csv";

    bool escribe_cabecera = false;
    if (!std::filesystem::exists(path) || std::filesystem::file_size(path) == 0) {
        escribe_cabecera = true;
    }

    std::ofstream csv_file;
    csv_file.open(path, std::ios::app); 

    if (!csv_file) {
        std::cerr << "Error: No se pudo abrir el archivo en resultados/res.csv" << std::endl;
        return 1;
    }

    if (escribe_cabecera) {
        csv_file << "nombre;steps;error_relativo;tiempo_ms\n";
    }

    // Definimos los pasos de subdivisión para todos los tests
    std::vector<std::size_t> steps = {32, 128, 512, 2048, 4096, 8192, 16384, 32768, 65636};
    const std::size_t Dim = 2; // Dimensiones de Integración

    std::cout << "Hardware concurrency: " 
              << std::thread::hardware_concurrency() 
              << " threads" << std::endl;

    // Parámetros comunes
    std::array<float, Dim> a; a.fill(5.0f);
    std::array<float, Dim> u; u.fill(0.5f);
    std::array<float, Dim> c; c.fill(5.0f);
    std::array<float, Dim> w_vec; w_vec.fill(0.5f);
    float w_scalar = 0.5f;

    // Tests
    run_test<GenzContinuous<Dim>, Dim>("Continuous", GenzContinuous<Dim>(a, u), steps, csv_file, run_name);
    run_test<GenzCornerpeak<Dim>, Dim>("Corner Peak", GenzCornerpeak<Dim>(c), steps, csv_file, run_name);
    run_test<GenzDiscontinuous<Dim>, Dim>("Discontinuous", GenzDiscontinuous<Dim>(a, u), steps, csv_file, run_name);
    run_test<GenzGaussianpeak<Dim>, Dim>("Gaussian Peak", GenzGaussianpeak<Dim>(w_vec, c), steps, csv_file, run_name);
    run_test<GenzOscilatory<Dim>, Dim>("Oscillatory", GenzOscilatory<Dim>(w_scalar, a), steps, csv_file, run_name);
    run_test<GenzProductpeak<Dim>, Dim>("Product Peak", GenzProductpeak<Dim>(w_vec, c), steps, csv_file, run_name);

    csv_file.close();
    return 0;
}