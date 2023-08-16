#pragma once
#include <functional>
#include <numbers>
#include <cmath>
#include <boost/math/special_functions/gamma.hpp>

namespace special_functions {

class GammaTables {

    static const int N_gamma_values = 172;

    const std::vector<double> gamma_int_values;
    const std::vector<double> gamma_half_int_values;

    std::vector<double> getGammaIntValues() const;
    std::vector<double> getGammaHalfIntValues() const;

    public:
        GammaTables();
        GammaTables(const GammaTables& old);
        double Gamma(int x) const;
        double GammaHalf(int x_minus_half) const;

};



// void storeFactorialValues();

// double getFactorial(int n);
// double getNegativeGammaHalf(int x_minus_half);



}