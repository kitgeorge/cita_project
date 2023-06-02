#pragma once
#include <functional>

namespace utility {

template<class Arg>
std::function<double(Arg)> 
addFunctions(std::function<double(Arg)> first,
             std::function<double(Arg)> second) {
    return [=] (Arg a) {
        return first(a) + second(a);
    };
}

template<class Arg0, class Arg1>
std::function<double(Arg0, Arg1)> 
addFunctions(std::function<double(Arg0, Arg1)> first,
             std::function<double(Arg0, Arg1)> second) {
    return [=] (Arg0 a0, Arg1 a1) {
        return first(a0, a1) + second(a0, a1);
    };
}

template<class Arg0, class Arg1, class Arg2>
std::function<double(Arg0, Arg1, Arg2)> 
addFunctions(std::function<double(Arg0, Arg1, Arg2)> first,
             std::function<double(Arg0, Arg1, Arg2)> second) {
    return [=] (Arg0 a0, Arg1 a1, Arg2 a2) {
        return first(a0, a1, a2) + second(a0, a1, a2);
    };
}

template <std::size_t N, class Arg0, class Arg1, class Arg2>
std::function<std::array<double, N>(Arg0, Arg1, Arg2)>
addFunctions(std::function<std::array<double, N>(Arg0, Arg1, Arg2)> first,
             std::function<std::array<double, N>(Arg0, Arg1, Arg2)> second) {
    return [=] (Arg0 a0, Arg1 a1, Arg2 a2) {
        std::array<double, N> first_array = first(a0, a1, a2);
        std::array<double, N> second_array = second(a0, a1, a2);
        std::array<double, N> output;
        for(int i = 0; i < N; ++i) {
            output[i] = first_array[i] + second_array[i];
        }
        return output;
    };
}

}