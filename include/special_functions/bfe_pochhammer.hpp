#pragma once
#include <vector>
#include <cmath>
#include <cassert>

namespace special_functions {
// Tabulates Pochhammer symbols in ranges required by BFE, given 
// constraints on n and l due to maximum calculable Gamma
double Pochhammer(double a, int i);

std::vector<std::vector<double>> getPochhammerIntValues();
std::vector<std::vector<double>> getPochhammerNegativeIntValues();
std::vector<std::vector<double>> getPochhammerHalfIntValues();
namespace {
// The names of these constants aren't quite precise,
// see how they're used instead
// Have to wrap the constants up in functions to avoid
// initialisation order trouble
int& int_a_max() {static int x = 100; return x;};
int& int_i_max() {static int x = 80; return x;};
int& negative_a_min() {static int x = -80; return x;};
int& negative_i_max() {static int x = 80; return x;};
int& half_a_max() {static int x = 171; return x;};
int& half_i_max() {static int x = 80; return x;};

// I'm not sure if there need to be pointers or something
// involved here
std::vector<std::vector<double>>& Pochhammer_int_values() {
    static std::vector<std::vector<double>>
    x = getPochhammerIntValues();
    return x;
}
std::vector<std::vector<double>>& Pochhammer_negative_int_values() {
    static std::vector<std::vector<double>>
    x = getPochhammerNegativeIntValues();
    return x;
}
std::vector<std::vector<double>>& Pochhammer_half_int_values() {
    static std::vector<std::vector<double>>
    x = getPochhammerHalfIntValues();
    return x;
}
}

double getPochhammerInt(int a, int i);
double getPochhammerHalfInt(int a_minus_half, int i);

}

