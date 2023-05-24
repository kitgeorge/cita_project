#include <vector>
#include <csignal>
#include <iostream>

#pragma once;

namespace utility {

template <typename T>
std::vector<std::vector<T>> makeShape(std::array<int, 2> shape) {
    std::vector<std::vector<T>>
    output(shape[0], std::vector<T>(shape[1]));
    return output;
}

template <typename T>
std::vector<std::vector<std::vector<T>>>
makeShape(std::array<int, 3> shape) {
    std::vector<std::vector<std::vector<T>>>
    output(shape[0], 
           std::vector<std::vector<T>>(shape[1],
                                       std::vector<T>(shape[2])));
    return output;
}

template <typename T>
std::vector<std::vector<std::vector<std::vector<T>>>>
makeShape(std::array<int, 4> shape) {
    std::vector<std::vector<std::vector<std::vector<T>>>>
    output(shape[0],
           std::vector<std::vector<std::vector<T>>>(shape[1],
                            std::vector<std::vector<T>>(shape[2], 
                                        std::vector<T>(shape[3]))));
    return output;
}

template <typename T>
std::vector<std::vector<std::vector<std::vector<T>>>>
setShape(std::vector<std::vector<std::vector<std::vector<T>>>> x,
         std::array<int, 4> shape) {
    x.resize(shape[0]);
    for(int i = 0; i < shape[0]; ++i) {
        x[i].resize(shape[1]);
        for(int j = 0; j < shape[1]; ++j) {
            x[i][j].resize(shape[2]);
            for(int k = 0; k < shape[2]; ++k) {
                x[i][j][k].resize(shape[3]);
            }
        }
    }
    return x;
}

template <typename T>
int getShape(std::vector<T> x) {
    return x.size();
}

template <typename T>
std::array<int, 2> getShape(std::vector<std::vector<T>> x) {
    std::array<int, 2> output;
    output[0] = x.size();
    output[1] = x[0].size();
    for(int i = 1; i < output[0]; ++i) {
        if(x[i].size() != output[1]) {
            std::cout << "Non-regular vector shape" << std::endl;
            std::abort();
        }
    }
    return output;
}

template <typename T>
std::array<int, 3> getShape(std::vector<std::vector<std::vector<T>>> x) {
    std::array<int, 3> output;
    output[0] = x.size();
    output[1] = x[0].size();
    output[2] = x[0][0].size();
    for(int i = 0; i < output[0]; ++i) {
        if(x[i].size() != output[1]) {
            std::cout << "Non-regular vector shape" << std::endl;
            std::abort();
        }
        for(int j = 0; j < output[1]; ++j) {
            if(x[i][j].size() != output[2]) {
                std::cout << "Non-regular vector shape" << std::endl;
                std::abort();
            }
        }
    }
    return output;
}

template <typename T>
std::array<int, 4> getShape(std::vector<std::vector<std::vector<std::vector<T>>>> x) {
    std::array<int, 4> output;
    output[0] = x.size();
    output[1] = x[0].size();
    output[2] = x[0][0].size();
    output[3] = x[0][0][0].size();
    for(int i = 0; i < output[0]; ++i) {
        if(x[i].size() != output[1]) {
            std::cout << "Non-regular vector shape" << std::endl;
            std::abort();
        }
        for(int j = 0; j < output[1]; ++j) {
            if(x[i][j].size() != output[2]) {
                std::cout << "Non-regular vector shape" << std::endl;
                std::abort();
            }
            for(int k = 0; k < output[2]; ++k) {
                if(x[i][j][k].size() != output[3]) {
                    std::cout << "Non-regualr vector shape" << std::endl;
                    std::abort();
                }
            }
        }
    }
    return output;
}


}