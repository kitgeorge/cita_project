#include <array>

namespace utility {

template <typename T, std::size_t N>
std::array<T, N> add_arrays(std::array<T, N> array1, std::array<T, N> array2) {
    std::array<T, N> output;
    for(int i = 0; i < N; ++i) {
        output[i] = array1[i] + array2[i];
    }
    return output;
}

template <typename T, std::size_t N_0, std::size_t N_1>
std::array<std::array<T, N_1>, N_0>
add_arrays(std::array<std::array<T, N_1>, N_0> array1,
           std::array<std::array<T, N_1>, N_0> array2) {
    std::array<std::array<T, N_1>, N_0> output;
    for(int i = 0; i < N_0; ++i) {
        for(int j = 0; j < N_1; ++j) {
            output[i][j] = array1[i][j] + array2[i][j];
        }
    }
    return output;
}

template <std::size_t N>
std::array<double, N> multiply_array(const std::array<double, N>& array, double scalar) {
    std::array<double, N> output;
    for(int i = 0; i < N; ++i) {
        output[i] = array[i]*scalar;
    }
    return output;
}

template <std::size_t N_0, std::size_t N_1>
std::array<std::array<double, N_1>, N_0> 
multiply_array(const std::array<std::array<double, N_1>, N_0>& array,
               double scalar) {
    std::array<std::array<double, N_1>, N_0>
    output;
    for(int i = 0; i < N_0; ++i) {
        for(int j = 0; j < N_1; ++j) {
            output[i][j] = array[i][j]*scalar;
        }
    }
    return output;
}
}