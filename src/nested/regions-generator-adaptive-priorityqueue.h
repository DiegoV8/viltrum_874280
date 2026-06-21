#pragma once
#include "nested.h"
#include "error-heuristic.h"
#include <mutex>
#include <optional>
#include <queue>
#include <vector>

// Adaptador para la Priority Queue con soporte thread-safe.
template <typename T, typename Comparator = std::less<T>>
class HeapAdapter {
private:
    std::priority_queue<T, std::vector<T>, Comparator> pq;
    mutable std::mutex mtx;

public:
    using value_type = T;

    void push(T val) {
        std::lock_guard<std::mutex> lock(mtx);
        pq.push(std::move(val));
    }

    void pop() {
        std::lock_guard<std::mutex> lock(mtx);
        if (!pq.empty())
            pq.pop();
    }

    T top() {
        std::lock_guard<std::mutex> lock(mtx);
        if (!pq.empty())
            return pq.top();
        return T();
    }

    bool empty() const {
        std::lock_guard<std::mutex> lock(mtx);
        return pq.empty();
    }

    std::optional<T> try_pop() {
        std::lock_guard<std::mutex> lock(mtx);
        if (!pq.empty()) {
            auto val = pq.top();
            pq.pop();
            return val;
        }
        return std::nullopt;
    }

    // Vuelca todo el contenido de la PQ a un vector y la deja vacía.
    std::vector<T> drain() {
        std::lock_guard<std::mutex> lock(mtx);
        std::vector<T> result;
        result.reserve(pq.size());
        while (!pq.empty()) {
            result.push_back(pq.top());
            pq.pop();
        }
        return result;
    }
};

// Comparador: ordena por error descendente (mayor error = mayor prioridad).
struct PqComparator {
    template<typename T>
    bool operator()(const T& a, const T& b) const {
        auto errA = std::get<0>(a.extra());
        auto errB = std::get<0>(b.extra());

        bool nanA = std::isnan(errA);
        bool nanB = std::isnan(errB);

        if (nanA && nanB) return false;
        if (nanA) return true;
        if (nanB) return false;

        return errA < errB; // max-heap: mayor error tiene prioridad
    }
};

namespace viltrum {

template<typename Rule, typename ErrorHeuristic,
         typename = std::enable_if_t<is_nested<Rule>::value>>
class RegionsGeneratorAdaptivePq {
    Rule rule;
    ErrorHeuristic error_heuristic;
    std::size_t subdivisions;

public:
    RegionsGeneratorAdaptivePq(const Rule& r, const ErrorHeuristic& er, std::size_t subdivisions)
        : rule(r), error_heuristic(er), subdivisions(subdivisions) {}

    template<std::size_t DIMBINS, typename F, typename Float, std::size_t DIM, typename Logger>
    auto generate(const std::array<std::size_t, DIMBINS>& bin_resolution,
                  const F& f, const Range<Float, DIM>& range, Logger& logger) const {

        auto r_init = region(f, rule, range.min(), range.max());
        auto errdim_init = error_heuristic(r_init);
        using ERegion = ExtendedRegion<decltype(r_init), decltype(errdim_init)>;

        HeapAdapter<ERegion, PqComparator> heap;
        heap.push(ERegion(r_init, errdim_init));

        std::atomic<std::size_t> completed{0};

        logger.log_progress(std::size_t(0), subdivisions);

        #pragma omp parallel
        {
            while (true) {
                std::size_t current = completed.load(std::memory_order_relaxed);
                if (current >= subdivisions) break;

                std::optional<ERegion> opt_r = heap.try_pop();
                if (!opt_r) continue;

                std::size_t my_slot = completed.fetch_add(1, std::memory_order_relaxed);
                if (my_slot >= subdivisions) {
                    heap.push(*opt_r);
                    break;
                }

                auto subregions = opt_r->split(f, std::get<1>(opt_r->extra()));
                for (auto& sr : subregions)
                    heap.push(ERegion(sr, error_heuristic(sr)));

                #pragma omp critical(logger_update)
                {
                    logger.log_progress(my_slot + 1, subdivisions);
                }
            }
        }

        logger.log_progress(subdivisions, subdivisions);

        return heap.drain();
    }
};

template<typename Rule, typename ErrorHeuristic>
auto regions_generator_adaptive_priorityqueue(const Rule& r, const ErrorHeuristic& er,
                                              std::size_t subdivisions) {
    return RegionsGeneratorAdaptivePq<Rule, ErrorHeuristic>(r, er, subdivisions);
}

} // namespace viltrum