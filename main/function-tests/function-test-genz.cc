#include <iostream>
#include <vector>
#include <chrono>
#include <iomanip>
#include <fstream>
#include <filesystem>

#include "../../viltrum.h"

using namespace viltrum;
using namespace viltrumtest;

/**
 * @brief Funcion para probar la integración de una función usando Priorityqueue.
 * @param name Nombre de la función a probar.
 * @param f Función a integrar.
 * @param steps Numero de pasos para la integración.
 * @param file Fichero csv de output.
 * @param run_name Nombre de la estructura usada.
 */
template<typename Function, std::size_t N>
void run_test_Pq(const std::string& name, const Function& f, const std::vector<std::size_t>& steps, std::ofstream& file, const std::string& run_name) {
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

/**
 * @brief Funcion para probar la integración de una función usando Multiqueue.
 * @param name Nombre de la función a probar.
 * @param f Función a integrar.
 * @param steps Numero de pasos para la integración.
 * @param file Fichero csv de output.
 * @param run_name Nombre de la estructura usada.
 */
template<typename Function, std::size_t N>
void run_test_Mq(const std::string& name, const Function& f, const std::vector<std::size_t>& steps, std::ofstream& file, const std::string& run_name) {
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
        auto integrator = integrator_adaptive_iterations_parallel_Mq(
            nested(simpson, trapezoidal), 
            s
        );

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

/**
 * @brief Funcion para probar la integración de una función usando Skiplist.
 * @param name Nombre de la función a probar.
 * @param f Función a integrar.
 * @param steps Numero de pasos para la integración.
 * @param file Fichero csv de output.
 * @param run_name Nombre de la estructura usada.
 */
template<typename Function, std::size_t N>
void run_test_Sl(const std::string& name, const Function& f, const std::vector<std::size_t>& steps, std::ofstream& file, const std::string& run_name) {
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
        auto integrator = integrator_adaptive_iterations_parallel_Sl(
            nested(simpson, trapezoidal), 
            s
        );

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
        std::cerr << "Uso: " << argv[0] << " <dataStructure>" << std::endl;
        return 1;
    }
    std::string run_name = "null";
    int ds = atoi(argv[1]);
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

    std::vector<std::size_t> steps = {32, 128, 512, 2048, 4096, 8192, 16384, 32768, 65636};
    const std::size_t Dim = 6;

    std::cout << "Hardware concurrency: " 
              << std::thread::hardware_concurrency() 
              << " threads" << std::endl;

    std::array<float, Dim> a; a.fill(5.0f);
    std::array<float, Dim> u; u.fill(0.5f);
    std::array<float, Dim> c; c.fill(5.0f);
    std::array<float, Dim> w_vec; w_vec.fill(0.5f);
    float w_scalar = 0.5f;

    switch (ds)
    {
    case 0:
        std::cout << "Running Priorityqueue tests..." << std::endl;
        run_name = "priorityqueue";
        run_test_Pq<GenzContinuous<Dim>, Dim>("Continuous", GenzContinuous<Dim>(a, u), steps, csv_file, run_name);
        run_test_Pq<GenzCornerpeak<Dim>, Dim>("Corner Peak", GenzCornerpeak<Dim>(c), steps, csv_file, run_name);
        run_test_Pq<GenzDiscontinuous<Dim>, Dim>("Discontinuous", GenzDiscontinuous<Dim>(a, u), steps, csv_file, run_name);
        run_test_Pq<GenzGaussianpeak<Dim>, Dim>("Gaussian Peak", GenzGaussianpeak<Dim>(w_vec, c), steps, csv_file, run_name);
        run_test_Pq<GenzOscilatory<Dim>, Dim>("Oscillatory", GenzOscilatory<Dim>(w_scalar, a), steps, csv_file, run_name);
        run_test_Pq<GenzProductpeak<Dim>, Dim>("Product Peak", GenzProductpeak<Dim>(w_vec, c), steps, csv_file, run_name);
        break;

    case 1:
        std::cout << "Running Multiqueue tests..." << std::endl;
        run_name = "multiqueue";
        run_test_Mq<GenzContinuous<Dim>, Dim>("Continuous", GenzContinuous<Dim>(a, u), steps, csv_file, run_name);
        run_test_Mq<GenzCornerpeak<Dim>, Dim>("Corner Peak", GenzCornerpeak<Dim>(c), steps, csv_file, run_name);
        run_test_Mq<GenzDiscontinuous<Dim>, Dim>("Discontinuous", GenzDiscontinuous<Dim>(a, u), steps, csv_file, run_name);
        run_test_Mq<GenzGaussianpeak<Dim>, Dim>("Gaussian Peak", GenzGaussianpeak<Dim>(w_vec, c), steps, csv_file, run_name);
        run_test_Mq<GenzOscilatory<Dim>, Dim>("Oscillatory", GenzOscilatory<Dim>(w_scalar, a), steps, csv_file, run_name);
        run_test_Mq<GenzProductpeak<Dim>, Dim>("Product Peak", GenzProductpeak<Dim>(w_vec, c), steps, csv_file, run_name);
        break;

    case 2:
        std::cout << "Running SkipList tests..." << std::endl;
        run_name = "skiplist";
        run_test_Sl<GenzContinuous<Dim>, Dim>("Continuous", GenzContinuous<Dim>(a, u), steps, csv_file, run_name);
        run_test_Sl<GenzCornerpeak<Dim>, Dim>("Corner Peak", GenzCornerpeak<Dim>(c), steps, csv_file, run_name);
        run_test_Sl<GenzDiscontinuous<Dim>, Dim>("Discontinuous", GenzDiscontinuous<Dim>(a, u), steps, csv_file, run_name);
        run_test_Sl<GenzGaussianpeak<Dim>, Dim>("Gaussian Peak", GenzGaussianpeak<Dim>(w_vec, c), steps, csv_file, run_name);
        run_test_Sl<GenzOscilatory<Dim>, Dim>("Oscillatory", GenzOscilatory<Dim>(w_scalar, a), steps, csv_file, run_name);
        run_test_Sl<GenzProductpeak<Dim>, Dim>("Product Peak", GenzProductpeak<Dim>(w_vec, c), steps, csv_file, run_name);
        break;

    case 3:
        std::cout << "Running all tests..." << std::endl;
        run_name = "multiqueue";
        std::cout << "Tests: " << run_name << std::endl;
        run_test_Mq<GenzContinuous<Dim>, Dim>("Continuous", GenzContinuous<Dim>(a, u), steps, csv_file, run_name);
        run_test_Mq<GenzCornerpeak<Dim>, Dim>("Corner Peak", GenzCornerpeak<Dim>(c), steps, csv_file, run_name);
        run_test_Mq<GenzDiscontinuous<Dim>, Dim>("Discontinuous", GenzDiscontinuous<Dim>(a, u), steps, csv_file, run_name);
        run_test_Mq<GenzGaussianpeak<Dim>, Dim>("Gaussian Peak", GenzGaussianpeak<Dim>(w_vec, c), steps, csv_file, run_name);
        run_test_Mq<GenzOscilatory<Dim>, Dim>("Oscillatory", GenzOscilatory<Dim>(w_scalar, a), steps, csv_file, run_name);
        run_test_Mq<GenzProductpeak<Dim>, Dim>("Product Peak", GenzProductpeak<Dim>(w_vec, c), steps, csv_file, run_name);
        run_name = "priorityqueue";
        std::cout << "Tests: " << run_name << std::endl;
        run_test_Pq<GenzContinuous<Dim>, Dim>("Continuous", GenzContinuous<Dim>(a, u), steps, csv_file, run_name);
        run_test_Pq<GenzCornerpeak<Dim>, Dim>("Corner Peak", GenzCornerpeak<Dim>(c), steps, csv_file, run_name);
        run_test_Pq<GenzDiscontinuous<Dim>, Dim>("Discontinuous", GenzDiscontinuous<Dim>(a, u), steps, csv_file, run_name);
        run_test_Pq<GenzGaussianpeak<Dim>, Dim>("Gaussian Peak", GenzGaussianpeak<Dim>(w_vec, c), steps, csv_file, run_name);
        run_test_Pq<GenzOscilatory<Dim>, Dim>("Oscillatory", GenzOscilatory<Dim>(w_scalar, a), steps, csv_file, run_name);
        run_test_Pq<GenzProductpeak<Dim>, Dim>("Product Peak", GenzProductpeak<Dim>(w_vec, c), steps, csv_file, run_name);
        run_name = "skiplist";
        std::cout << "Tests: " << run_name << std::endl;
        run_test_Sl<GenzContinuous<Dim>, Dim>("Continuous", GenzContinuous<Dim>(a, u), steps, csv_file, run_name);
        run_test_Sl<GenzCornerpeak<Dim>, Dim>("Corner Peak", GenzCornerpeak<Dim>(c), steps, csv_file, run_name);
        run_test_Sl<GenzDiscontinuous<Dim>, Dim>("Discontinuous", GenzDiscontinuous<Dim>(a, u), steps, csv_file, run_name);
        run_test_Sl<GenzGaussianpeak<Dim>, Dim>("Gaussian Peak", GenzGaussianpeak<Dim>(w_vec, c), steps, csv_file, run_name);
        run_test_Sl<GenzOscilatory<Dim>, Dim>("Oscillatory", GenzOscilatory<Dim>(w_scalar, a), steps, csv_file, run_name);
        run_test_Sl<GenzProductpeak<Dim>, Dim>("Product Peak", GenzProductpeak<Dim>(w_vec, c), steps, csv_file, run_name);
        break;
    
    default:
        std::cout << "Estructura de datos no reconocida" << std::endl;
        break;
    }

    csv_file.close();
    return 0;
}