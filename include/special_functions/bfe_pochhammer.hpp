#pragma once
#include <vector>
#include <cmath>
#include <cassert>
#include <boost/multiprecision/gmp.hpp>

using LooongDouble = boost::multiprecision::mpf_float_1000;

namespace special_functions {
// Tabulates Pochhammer symbols in ranges required by BFE, given 
// constraints on n and l due to maximum calculable Gamma
LooongDouble Pochhammer(double a, int i);

std::vector<std::vector<LooongDouble>> getPochhammerIntValues();
std::vector<std::vector<LooongDouble>> getPochhammerNegativeIntValues();
std::vector<std::vector<LooongDouble>> getPochhammerHalfIntValues();
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
std::vector<std::vector<LooongDouble>>& Pochhammer_int_values() {
    static std::vector<std::vector<LooongDouble>> x;
    // Not sure if the compiler would optimise for the calculation to
    // only happen once anyway, but this makes sure
    if(x.empty()) {x = getPochhammerIntValues();};
    return x;
}
std::vector<std::vector<LooongDouble>>& Pochhammer_negative_int_values() {
    static std::vector<std::vector<LooongDouble>> x;
    if(x.empty()) {x = getPochhammerNegativeIntValues();};
    return x;
}
std::vector<std::vector<LooongDouble>>& Pochhammer_half_int_values() {
    static std::vector<std::vector<LooongDouble>> x;
    if(x.empty()) {x = getPochhammerHalfIntValues();};
    return x;
}
}

LooongDouble getPochhammerInt(int a, int i);
LooongDouble getPochhammerHalfInt(int a_minus_half, int i);

}

