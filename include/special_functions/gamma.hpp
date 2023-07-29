#pragma once
#include <functional>
#include <numbers>
#include <cmath>
#include <boost/math/special_functions/gamma.hpp>

namespace special_functions {

namespace {
static const int N_factorial_values = 85;
static std::vector<double> factorial_values;
static std::vector<double> negative_gamma_values;
}

void storeFactorialValues();

double getFactorial(int n);
double getNegativeGammaHalf(int x_minus_half);


double Gamma(int x);
double GammaHalf(int x_minus_half);
double Pochhammer(double a, int i);


}