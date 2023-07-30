#include "gamma.hpp"
#include <cassert>
#include <cmath>
#include <iostream>

namespace special_functions {

std::vector<double> getGammaIntValues() {
    std::vector<double> output(N_gamma_values);
    // 0 undefined
    for(int i = 1; i < N_gamma_values; ++i) {
        output[i] = boost::math::tgamma(i);
    }
    return output;
}

std::vector<double> getGammaHalfIntValues() {
    std::vector<double> output(N_gamma_values);
    for(int i = 0; i < N_gamma_values; ++i) {
        output[i] = boost::math::tgamma(i + 0.5);
    }
    return output;
}

double Gamma(int x) {
    assert(x > 0);
    assert(x < N_gamma_values);
    return gamma_int_values[x];
}

double GammaHalf(int x) {
    assert(x >= 0);
    assert(x < N_gamma_values);
    return gamma_half_int_values[x];
}

// void storeFactorialValues() {
//     std::function<double(int)>
//     factorial = [] (int n) {
//         assert(n >= 0);
//         double output = 1;
//         for(int i = 2; i <= n; ++i) {
//             output = output*i;
//         }
//         assert(output > 0);
//         return output;
//     };
//     factorial_values.resize(N_factorial_values);
//     for(int i = 0; i < N_factorial_values; ++i) {
//         factorial_values[i] = factorial(i);
//     }
//     negative_gamma_values.resize(N_factorial_values);
//     for(int i = 0; i < N_factorial_values; ++i) {
//         negative_gamma_values[i] = boost::math::tgamma(-i - 0.5);
//     }
// }

// double getFactorial(int n) {
//     assert(n >= 0);
//     assert(n < N_factorial_values);
//     if(factorial_values.empty()) {
//         storeFactorialValues();
//     }
//     return factorial_values[n];
// }

// double getNegativeGammaHalf(int x_minus_half) {
//     assert(x_minus_half < 0);
//     assert(x_minus_half >= -N_factorial_values);
//     if(negative_gamma_values.empty()) {
//         storeFactorialValues();
//     }
//     return negative_gamma_values[-x_minus_half - 1];
// }


// double Gamma(int x) {
//     return getFactorial(x - 1);
// }

// double GammaHalf(int x_minus_half) {
//     if(x_minus_half < 0) {
//         std::cout << "NEGATIVE" << std::endl;
//         return getNegativeGammaHalf(x_minus_half);
//     }

//     double output = getFactorial(2*x_minus_half)/(pow(4, x_minus_half)
//                                       *getFactorial(x_minus_half))
//                     *sqrt(std::numbers::pi);
//     if(!std::isfinite(output)) {
//         std::cout << "Gamma error , " << x_minus_half << ", " << output << std::endl;
//     }
//     assert(output > 0);
//     assert(std::isfinite(output));
//     return output;
// }

double Pochhammer(double a, int i) {
    assert(i >= 0);
    double output = 1; 
    for(int j = 0; j < i; ++j) {
        output = output*(a + j);
    }
    assert(std::isfinite(output));
    return output;
}
}