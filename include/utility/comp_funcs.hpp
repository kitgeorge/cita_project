#pragma once
#include <functional>
#include <complex>
#include <array>

namespace utility {

template<class Arg0, class Arg1>
std::function<std::complex<double>(Arg0, Arg1)>
conjugateFunction(std::function<std::complex<double>(Arg0, Arg1)> func) {
    return [=] (Arg0 a0, Arg1 a1) {
        return std::conj(func(a0, a1));
    };
}

template<class Arg0, class Arg1, std::size_t N>
std::function<std::array<std::complex<double>, N>(Arg0, Arg1)>
conjugateFunction(std::function<std::array<std::complex<double>, N>(Arg0, Arg1)> func) {
    return [=] (Arg0 a0, Arg1 a1) {
        std::array<std::complex<double>, N> output;
        std::array<std::complex<double>, N> value = func(a0, a1);
        for(int i = 0; i < N; ++i) {
            output[i] = std::conj(value[i]);
        }
        return output;
    };
}

template<class Arg0, class Arg1>
std::function<double(Arg0, Arg1)>
realFunction(std::function<std::complex<double>(Arg0, Arg1)> func) {
    return [=] (Arg0 a0, Arg1 a1) {
        return func(a0, a1).real();
    };
}

template<class Arg0, class Arg1, std::size_t N>
std::function<std::array<double, N>(Arg0, Arg1)>
realFunction(std::function<std::array<std::complex<double>, N>(Arg0, Arg1)> func) {
    return [=] (Arg0 a0, Arg1 a1) {
        std::array<double, N> output;
        std::array<std::complex<double>, N> value = func(a0, a1);
        for(int i = 0; i < N; ++i) {
            output[i] = value[i].real();
        }
        return output;
    };
}


}