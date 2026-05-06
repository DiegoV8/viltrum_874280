#pragma once

//#include "function-Nd.h"

#include <array>
#include <cmath>
#include <tuple>

namespace viltrumtest {

// Taken from https://www.sfu.ca/~ssurjano/integration.html
//   Analyitical solution here. https://www.sciencedirect.com/science/article/pii/S0893965924004646
//   I found a different one here: https://proceedings.neurips.cc/paper/2020/file/3fe94a002317b5f9259f82690aeea4cd-Supplemental.pdf
//   and here: https://github.com/ImperialCollegeLondon/BART-Int/blob/master/src/genz/analyticalIntegrals.R
template<std::size_t N>
class GenzCornerpeak /*: public pattern::Reflectable<GenzCornerpeak<N>,FunctionNDBase<N>>*/ {
    std::array<float,N> c;
public:
    GenzCornerpeak(const std::array<float,N>& c) : c(c) {}
    GenzCornerpeak() { c.fill(5.0f); } 
    static const char* type_name() { return "genz-cornerpeak"; }
    auto reflect() { return std::tie(c); }
    auto reflect_names() const { return std::tuple("c"); }
    float operator()(const std::array<float,N>& x) const /*override*/ {
        float sum=1.0f;
        for (std::size_t i = 0; i<N; ++i) sum += c[i]*x[i];
        return 1.0f/std::pow(sum,N+1);
    }
    
    /* This didn't work so we impemented the neurips one 
    float integral_primary() const override {
        float fact_d = 1.0f; for (std::size_t i = 1; i<=N; ++i) fact_d*=i;
        float prod_c = 1.0f; for (std::size_t i = 0; i<N; ++i) prod_c*=c[i];
        float sum = 1.0;
        for (std::size_t i = 0; i < N; ++i) {
            float sum_local = 0.0f;
            for (std::size_t j = 0; j<=i; ++j) {
                sum_local += c[j]; 
            }
            sum += std::pow(-1.0f,i+1)/(1.0f + sum_local); 
        } 
        return sum/(fact_d*prod_c);
    } 
    */
private:
    static float factorial(unsigned long n) {
        float f(1);
        for (unsigned long i=2;i<=n;++i) f*=float(i);
        return f;
    }
    static float choose(unsigned long n, unsigned long k) {
        return factorial(n)/(factorial(k)*factorial(n-k));
    }  
public:
    //This below assumes that all values are equal but the previous one didn't work...
    float integral_primary() const /*override*/ {
        float prod_c = 1.0f; for (std::size_t i = 0; i<N; ++i) prod_c*=c[i];
        float sum = 0.0;
        for (std::size_t i = 1; i <= N; ++i) 
            sum += choose(N-1,i-1)*std::pow(-1.0f,i-1)/((c[i-1]*i + 1.0f)*(c[i-1]*(i-1) + 1.0f));
        
        return c[0]*sum/(factorial(N)*prod_c); //c[0] is invented because we assume all of them are equal 
    } 

};



}