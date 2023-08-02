#include "bfe_pochhammer.hpp"
#include <iostream>

namespace special_functions {


LooongDouble Pochhammer(double a, int i) {
    assert(i >= 0);
    assert(std::isfinite(a));
    assert(std::isfinite(i));
    LooongDouble output = 1; 
    for(int j = 0; j < i; ++j) {
        output = output*(a + j);
    }
    // std::cout << a << ", " << i << ", " << output << std::endl;
    assert(std::isfinite(output));
    return output;
}

std::vector<std::vector<double>> getPochhammerIntValues() {
    std::vector<std::vector<double>>
    output(int_a_max() + 1, std::vector<double>(int_i_max() + 1));
    for(int i = 0; i <= int_a_max(); ++i) {
        for(int j = 0; j <= int_i_max(); ++j) {
            output[i][j] = Pochhammer(i, j);
        }
    }
    return output;
}

std::vector<std::vector<LooongDouble>> getPochhammerNegativeIntValues() {
    std::vector<std::vector<LooongDouble>>
    output(-negative_a_min() + 1, std::vector<LooongDouble>(negative_i_max() + 1));
    for(int i = 0; i >= negative_a_min(); --i) {
        for(int j = 0; j <= negative_i_max(); ++j) {
            output[-i][j] = Pochhammer(i, j);
        }
    }
    return output;
}

std::vector<std::vector<LooongDouble>> getPochhammerHalfIntValues() {
    std::vector<std::vector<LooongDouble>>
    output(half_a_max() + 1, std::vector<LooongDouble>(half_i_max() + 1));
    for(int i = 0; i <= half_a_max(); ++i) {
        for(int j = 0; j <= half_i_max(); ++j) {
            output[i][j] = Pochhammer(i + 0.5, j);
        }
    }
    return output;
}


LooongDouble getPochhammerInt(int a, int i) {
    assert(i >= 0);
    if(a >= 0) {
        assert(a <= int_a_max());
        assert(i <= int_i_max());
        // std::cout << a << ", " << i << ", " << Pochhammer_int_values()[a][i] << std::endl;
        return Pochhammer_int_values()[a][i];
    }
    else {
        assert(a >= negative_a_min());
        assert(i <= negative_i_max());
        // std::cout << a << ", " << i << ", " << Pochhammer_negative_int_values()[-a][i] << std::endl;
        return Pochhammer_negative_int_values()[-a][i];
    }
}

LooongDouble getPochhammerHalfInt(int a_minus_half, int i) {
    assert(a_minus_half >= 0);
    assert(a_minus_half <= half_a_max());
    assert(i >= 0);
    assert(i <= half_i_max());
    // std::cout << "Half, " << a_minus_half << ", " << i << ", " << Pochhammer_half_int_values()[a_minus_half][i] << std::endl;
    return Pochhammer_half_int_values()[a_minus_half][i];
}

}