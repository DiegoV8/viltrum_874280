#include "../../viltrum.h"
#include <iostream>
#include <iomanip>

int main() {

    float sol = viltrum::integrate(
        viltrum::rectangle_rule(64),
        [] (float x) -> float { return std::sin(x);}, 
        viltrum::range(0.0f, 3.14159265f)             
    );

    std::cout<<std::fixed<<std::setprecision(6)<<sol<<" should be close to 2\n";
}