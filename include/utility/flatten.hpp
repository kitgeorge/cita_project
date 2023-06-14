#include <vector>
#include <iostream>
#include <cstdlib>
#include "shape.hpp"

namespace utility {

template <typename T>
std::vector<T> flatten(std::vector<std::vector<T>> data) {
    std::array<int, 2> dimensions = getShape(data);

    std::vector<T> output(dimensions[0]*dimensions[1]);
    for(int i = 0; i < dimensions[0]; ++i) {
        for(int j = 0; j < dimensions[1]; ++j) {
            output[i*dimensions[1] + j] = data[i][j];
        }
    }
    return output;
}

template <typename T>
std::vector<T> flatten(std::vector<std::vector<std::vector<T>>> data) {
    std::array<int, 3> dimensions = getShape(data);
 
    std::vector<T> output(dimensions[0]*dimensions[1]*dimensions[2]);
    for(int i = 0; i < dimensions[0]; ++i) {
        for(int j = 0; j < dimensions[1]; ++j) {
            for(int k = 0; k < dimensions[2]; ++k) {
                output[i*dimensions[1]*dimensions[2]
                       + j*dimensions[2] + k] = data[i][j][k];
            }
        }
    }
    return output;
}

template <typename T>
std::vector<T> flatten(std::vector<std::vector<std::vector<std::vector<T>>>> data) {
    std::array<int, 4> dimensions = getShape(data);

    std::vector<T> output(dimensions[0]*dimensions[1]*dimensions[2]*dimensions[3]);
    for(int i = 0; i < dimensions[0]; ++i) {
        for(int j = 0; j < dimensions[1]; ++j) {
            for(int k = 0; k < dimensions[2]; ++k) {
                for(int l = 0; l < dimensions[3]; ++l) {
                    output[i*dimensions[1]*dimensions[2]*dimensions[3]
                        + j*dimensions[2]*dimensions[3] 
                        + k*dimensions[3] + l] = data[i][j][k][l];
                    
                }
            }
        }
    }
    return output;
}

template <typename T, std::size_t N>
std::vector<T> flatten(std::array<T, N> data) {
    std::vector<T> output(N);
    for(int i = 0; i < N; ++i) {
        output[i] = data[i];
    }
    return output;
}

template <typename T, std::size_t N>
std::vector<T> flatten(std::vector<std::array<T, N>> data) {
    int N_0 = data.size();
    std::vector<T> output(N_0*N);
    for(int i = 0; i < N_0; ++i) {
        for(int j = 0; j < N; ++j) {
            output[i*N + j] = data[i][j];
        }
    }
    return output;
}


}