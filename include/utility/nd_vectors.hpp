#pragma once
#include <vector>
#include <array>
#include <cassert>

namespace utility {

template <typename T>
using vector2d = std::vector<std::vector<T>>;

template <typename T>
using vector3d = std::vector<std::vector<std::vector<T>>>;

template <typename T>
using vector4d = std::vector<std::vector<std::vector<std::vector<T>>>>;

template <typename T>
using vector5d = std::vector<std::vector<std::vector<std::vector<std::vector<T>>>>>;


template <typename T>
vector2d<T> makeShape(std::array<int, 2> shape) {
    vector2d<T> output(shape[0], std::vector<T>(shape[1]));
    return output;
}

template <typename T>
vector3d<T> makeShape(std::array<int, 3> shape) {
    std::array<int, 2> sub_shape = {{shape[1], shape[2]}};
    vector3d<T> output(shape[0], makeShape<T>(sub_shape));
    return output;
}

template <typename T>
vector4d<T> makeShape(std::array<int, 4> shape) {
    std::array<int, 3> sub_shape = {{shape[1], shape[2], shape[3]}};
    vector4d<T> output(shape[0], makeShape<T>(sub_shape));
    return output;
}

template <typename T>
vector5d<T> makeShape(std::array<int, 5> shape) {
    std::array<int, 4> sub_shape = {{shape[1], shape[2], shape[3], shape[4]}};
    vector5d<T> output(shape[0], makeShape<T>(sub_shape));
    return output;
}


template <typename T>
vector2d<T> reshape(std::vector<T> flat, std::array<int, 2> shape) {
    assert(flat.size() == shape[0]*shape[1]);
    vector2d<T> output = makeShape<T>(shape);
    for(int i = 0; i < shape[0]; ++i) {
        for(int j = 0; j < shape[1]; ++j) {
            output[i][j] = flat[i*shape[1] + j];
        }
    }
    return output;
}

template <typename T>
vector3d<T> reshape(std::vector<T> flat, std::array<int, 3> shape) {
    assert(flat.size() == shape[0]*shape[1]*shape[2]);
    vector3d<T> output = makeShape<T>(shape);
    for(int i = 0; i < shape[0]; ++i) {
        for(int j = 0; j < shape[1]; ++j) {
            for(int k = 0; k < shape[2]; ++k) {
                output[i][j][k] = flat[i*shape[1]*shape[2] + j*shape[2] + k];
            }
        }
    }
    return output;
}

template <typename T>
vector4d<T> reshape(std::vector<T> flat, std::array<int, 4> shape) {
    assert(flat.size() == shape[0]*shape[1]*shape[2]*shape[3]);
    vector4d<T> output = makeShape<T>(shape);
    for(int i = 0; i < shape[0]; ++i) {
        for(int j = 0; j < shape[1]; ++j) {
            for(int k = 0; k < shape[2]; ++k) {
                for(int l = 0; l < shape[3]; ++l) {
                    output[i][j][k][l] = flat[i*shape[1]*shape[2]*shape[3] 
                                              + j*shape[2]*shape[3] 
                                              + k*shape[3] + l];
                }
            }
        }
    }
    return output;
}


}