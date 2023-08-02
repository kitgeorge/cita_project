#pragma once
#include <functional>
#include <numbers>
#include <cmath>
#include <boost/math/special_functions/gamma.hpp>

namespace special_functions {

std::vector<double> getGammaIntValues();
std::vector<double> getGammaHalfIntValues();

namespace {
// static const int N_factorial_values = 85;
// static std::vector<double> factorial_values;
// static std::vector<double> negative_gamma_values;

// Tabulates Gamma function for all integers and half-integers 0 < x < 172
// (at compile time!)
int& N_gamma_values() {static int x = 172; return x;}

static std::vector<double>& gamma_int_values() {
    static std::vector<double> x;
    if(x.empty()) {x = getGammaIntValues();};
    return x;
}
static std::vector<double>& gamma_half_int_values() {
    static std::vector<double> x;
    if(x.empty()) {x = getGammaHalfIntValues();};
    return x;
}
// static const std::vector<double> gamma_int_values = getGammaIntValues();
// static const std::vector<double> gamma_half_int_values = getGammaHalfIntValues();
}

// void storeFactorialValues();

// double getFactorial(int n);
// double getNegativeGammaHalf(int x_minus_half);


double Gamma(int x);
double GammaHalf(int x_minus_half);


}