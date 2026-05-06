#pragma once

//#include "function-Nd.h"

#include <array>
#include <cmath>
#include <tuple>

namespace viltrumtest {

// Taken from https://www.sfu.ca/~ssurjano/integration.html
//   Analyitical solution here. https://www.sciencedirect.com/science/article/pii/S0893965924004646
template<std::size_t N>
class GenzOscilatory /*: public pattern::Reflectable<GenzOscilatory<N>,FunctionNDBase<N>>*/ {
    float w;
    std::array<float,N> c;
public:
    GenzOscilatory(float w, const std::array<float,N>& c) : w(w), c(c) {}
    GenzOscilatory() : w(0.5f) { c.fill(5.0f); } 
    static const char* type_name() { return "genz-oscilatory"; }
    auto reflect() { return std::tie(w,c); }
    auto reflect_names() const { return std::tuple("w","c"); }
    float operator()(const std::array<float,N>& x) const /*override*/ {
        float sum=0.0f;
        for (std::size_t i = 0; i<N; ++i) sum += c[i]*x[i];
        return std::cos(2*M_PI*w + sum);   
    }
    float integral_primary() const /*override*/ {
        float sum = 0.0;  for (std::size_t i = 0; i<N; ++i) sum += c[i];
        float prod = 1.0; for (std::size_t i = 0; i<N; ++i) prod *= ((std::sin(c[i]/2.0))/c[i]);
        return std::pow(2.0,N)*std::cos(2*M_PI*w + 0.5*sum)*prod; 
    } 
};



}