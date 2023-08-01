#include <vector>
#include <array>

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

}