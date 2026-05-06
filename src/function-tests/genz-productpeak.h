#pragma once

//#include "function-Nd.h"

#include <array>
#include <cmath>
#include <tuple>

namespace viltrumtest {

// Taken from https://www.sfu.ca/~ssurjano/integration.html
//   Analyitical solution here. https://www.sciencedirect.com/science/article/pii/S0893965924004646
template<std::size_t N>
class GenzProductpeak /*: public pattern::Reflectable<GenzProductpeak<N>,FunctionNDBase<N>>*/ {
    std::array<float,N> w, c;
public:
    GenzProductpeak(const std::array<float,N>& w, const std::array<float,N>& c) : w(w), c(c) {}
    GenzProductpeak() { w.fill(0.5), c.fill(5.0f); } 
    static const char* type_name() { return "genz-productpeak"; }
    auto reflect() { return std::tie(w,c); }
    auto reflect_names() const { return std::tuple("w","c"); }
    float operator()(const std::array<float,N>& x) const /*override*/ {
        float prod=1.0f;
        for (std::size_t i = 0; i<N; ++i) prod /= ((1.0f/(c[i]*c[i])) + (x[i] - w[i])*(x[i]-w[i]));
        return prod;
    }
    float integral_primary() const /*override*/ {
        float prod = 1.0f; 
        for (std::size_t i = 0; i<N; ++i) prod*=(c[i]*(std::atan(c[i]*w[i]) + std::atan(c[i] - c[i]*w[i])));
        return prod;
    } 
};



}