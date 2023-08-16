#pragma once
#include <vector>
#include <cmath>
#include <cassert>
#include <boost/multiprecision/gmp.hpp>

using LooongDouble = boost::multiprecision::mpf_float_1000;

namespace special_functions {
// Tabulates Pochhammer symbols in ranges required by BFE, given 
// constraints on n and l due to maximum calculable Gamma

class PochhammerTables {
    static const int int_a_max = 100;
    static const int int_i_max = 80;
    static const int negative_a_min = -80;
    static const int negative_i_max = 80;
    static const int half_a_max = 171;
    static const int half_i_max = 80;

    const std::vector<std::vector<LooongDouble>> 
    Pochhammer_int_values;
    const std::vector<std::vector<LooongDouble>> 
    Pochhammer_negative_int_values;
    const std::vector<std::vector<LooongDouble>> 
    Pochhammer_half_int_values;

    std::vector<std::vector<LooongDouble>> 
    getPochhammerIntValues() const;
    std::vector<std::vector<LooongDouble>> 
    getPochhammerNegativeIntValues() const;
    std::vector<std::vector<LooongDouble>> 
    getPochhammerHalfIntValues() const;

    public:
        PochhammerTables();
        PochhammerTables(const PochhammerTables& old);
        LooongDouble getPochhammerInt(int a, int i) const;
        LooongDouble getPochhammerHalfInt(int a_minus_half, int i) const;
};

}

