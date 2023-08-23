#pragma once
#include <functional>
#include <iostream>
#include <complex>

namespace utility {

// template<class Arg>
// std::function<double(Arg)> 
// addFunctions(std::function<double(Arg)> first,
//              std::function<double(Arg)> second) {
//     return [=] (Arg a) {
//         return first(a) + second(a);
//     };
// }

// template<class Arg0, class Arg1>
// std::function<double(Arg0, Arg1)> 
// addFunctions(std::function<double(Arg0, Arg1)> first,
//              std::function<double(Arg0, Arg1)> second) {
//     return [=] (Arg0 a0, Arg1 a1) {
//         return first(a0, a1) + second(a0, a1);
//     };
// }

// template<class Arg0, class Arg1, class Arg2>
// std::function<double(Arg0, Arg1, Arg2)> 
// addFunctions(std::function<double(Arg0, Arg1, Arg2)> first,
//              std::function<double(Arg0, Arg1, Arg2)> second) {
//     return [=] (Arg0 a0, Arg1 a1, Arg2 a2) {
//         return first(a0, a1, a2) + second(a0, a1, a2);
//     };
// }

// template <std::size_t N, class Arg0, class Arg1, class Arg2>
// std::function<std::array<double, N>(Arg0, Arg1, Arg2)>
// addFunctions(std::function<std::array<double, N>(Arg0, Arg1, Arg2)> first,
//              std::function<std::array<double, N>(Arg0, Arg1, Arg2)> second) {
//     return [=] (Arg0 a0, Arg1 a1, Arg2 a2) {
//         std::array<double, N> first_array = first(a0, a1, a2);
//         std::array<double, N> second_array = second(a0, a1, a2);
//         std::array<double, N> output;
//         for(int i = 0; i < N; ++i) {
//             output[i] = first_array[i] + second_array[i];
//         }
//         return output;
//     };
// }


// This is a mess, should sort out sometime


double addFunctions(const std::vector<std::function<double(double, double, double)>>& 
                    functions, double a0, double a1, double a2);

template <class Arg0, class Arg1>
std::function<std::complex<double>(Arg0, Arg1)>
addFunctions(const std::vector<std::function<std::complex<double>(Arg0, Arg1)>>&
             functions) {
    return [=] (Arg0 a0, Arg1 a1) {
        std::complex<double> output = 0;
        for(auto func: functions) {
            output += func(a0, a1);
        }
        return output;
    };
}

template <class Arg0, class Arg1, std::size_t N>
std::function<std::array<std::complex<double>, N>(Arg0, Arg1)>
addFunctions(const std::vector<std::function<std::array<std::complex<double>, N>(Arg0, Arg1)>>&
             functions) {
    return [=] (Arg0 a0, Arg1 a1) {
        std::array<std::complex<double>, N> output = {};
        for(auto func: functions) {
            std::array<std::complex<double>, N> value = func(a0, a1);
            for(int i = 0; i < N; ++i) {
                output[i] += value[i];
            }
        }
        return output;
    };
}

template <std::size_t N>
std::array<double, N> 
addFunctions(const std::vector<std::function<std::array<double, N>(double, double, double)>>& 
             functions, double a0, double a1, double a2) {
    std::array<double, N> output;
    for(int i = 0; i < N; ++i) {
        output[i] = 0;
    }
    for(auto f: functions) {
        std::array<double, N> temp = f(a0, a1, a2);
        for(int i = 0; i < N; ++i) {
            output[i] += temp[i];
        }
    }
    return output;
}


template<class Arg>
std::function<double(Arg)>
multiplyFunction(std::function<double(Arg)> f, double scalar) {
    return [=] (Arg a) {
        return scalar*f(a);
    };
}

template<class Arg0, class Arg1>
std::function<double(Arg0, Arg1)>
multiplyFunction(std::function<double(Arg0, Arg1)> f, double scalar) {
    return [=] (Arg0 a0, Arg1 a1) {
        return scalar*f(a0, a1);
    };
}

template<class Arg0, class Arg1>
std::function<std::complex<double>(Arg0, Arg1)>
multiplyFunction(const std::function<std::complex<double>(Arg0, Arg1)>& f, 
                 std::complex<double> scalar) {
    return [=] (Arg0 a0, Arg1 a1) {
        return scalar*f(a0, a1);
    };
}

template<class Arg0, class Arg1, std::size_t N>
std::function<std::array<std::complex<double>, N>(Arg0, Arg1)>
multiplyFunction(const std::function<std::array<std::complex<double>, N>(Arg0, Arg1)>& f, 
                 std::complex<double> scalar) {
    return [=] (Arg0 a0, Arg1 a1) {
        std::array<std::complex<double>, N> output = f(a0, a1);
        for(int i = 0; i < N; ++i) {
            output[i] = scalar*output[i];
        }
        return output;
    };
}

template<class Arg0, class Arg1, class Arg2>
std::function<double(Arg0, Arg1, Arg2)>
multiplyFunction(std::function<double(Arg0, Arg1, Arg2)> f, double scalar) {
    return [=] (Arg0 a0, Arg1 a1, Arg2 a2) {
        return scalar*f(a0, a1, a2);
    };
}

template <std::size_t N>
std::function<std::array<double, N>(double, double, double)>
multiplyFunction(const std::function<std::array<double, N>(double, double, double)>& 
             function, double scalar) {
    return [=] (double a0, double a1, double a2) {
        return scalar*function(a0, a1, a2);
    };
}