#pragma once

//#include "function-Nd.h"

#include <array>
#include <cmath>
#include <tuple>

namespace viltrumtest {

// Taken from https://www.sfu.ca/~ssurjano/integration.html
//   Analyitical solution here. https://www.sciencedirect.com/science/article/pii/S0893965924004646
//   or better from here: https://proceedings.neurips.cc/paper/2020/file/3fe94a002317b5f9259f82690aeea4cd-Supplemental.pdf 
template<std::size_t N>
class GenzDiscontinuous /*: public pattern::Reflectable<GenzDiscontinuous<N>,FunctionNDBase<N>>*/ {
    std::array<float,N> a;
    std::array<float,N> u;
public:
    GenzDiscontinuous(const std::array<float,N>& a, const std::array<float,N>& u) : a(a), u(u) {}
    GenzDiscontinuous(){ a.fill(5.0f); u.fill(0.5f); } 
    static const char* type_name() { return "genz-discontinuous"; }
    auto reflect() { return std::tie(a,u); }
    auto reflect_names() const { return std::tuple("a","u"); }
    float operator()(const std::array<float,N>& x) const /*override*/ {
        float sum=0.0f;
        for (std::size_t i = 0; i<N; ++i) {
            if (x[i]>u[i]) return 0; 
            sum += a[i]*x[i]; 
        } 
        return std::exp(sum);
    }
    
    float integral_primary() const /*override*/ {
        if constexpr(N==1) {
            return (std::exp(a[0]*u[0]) - 1.0f)/a[0];
        } else {
            float sol(1);
            for (std::size_t i = 0;i<N;++i) {
                sol *= (std::exp(a[i]*std::min(1.0f,u[i]))-1.0f)/a[i]; 
            } 
            return sol;
        } 
    } 
    
};



}