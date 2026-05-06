#pragma once

//#include "function-Nd.h"

#include <array>
#include <cmath>
#include <tuple>

namespace viltrumtest {

// Taken from https://www.sfu.ca/~ssurjano/integration.html
//   Analyitical solution here. https://www.sciencedirect.com/science/article/pii/S0893965924004646
template<std::size_t N>
class GenzGaussianpeak /*: public pattern::Reflectable<GenzGaussianpeak<N>,FunctionNDBase<N>>*/ {
    std::array<float,N> w, c;
public:
    GenzGaussianpeak(const std::array<float,N>& w, const std::array<float,N>& c) : w(w), c(c) {}
    GenzGaussianpeak() { w.fill(0.5), c.fill(5.0f); } 
    static const char* type_name() { return "genz-gaussianpeak"; }
    auto reflect() { return std::tie(w,c); }
    auto reflect_names() const { return std::tuple("w","c"); }
    float operator()(const std::array<float,N>& x) const /*override*/ {
        float sum=0.0f;
        for (std::size_t i = 0; i<N; ++i) sum += c[i]*c[i]*(x[i]-w[i])*(x[i]-w[i]);
        return std::exp(-sum);
    }
    float integral_primary() const /*override*/ {
        float prod = 1.0f; 
        for (std::size_t i = 0; i<N; ++i) prod*=(std::erf(c[i]*w[i]) + std::erf(c[i] - c[i]*w[i]))/c[i];
        return std::pow(M_PI,float(N)/2.0f)*prod/std::pow(2,N);
    } 
};



}