#include <functional>

namespace utility {

template<class R, class Arg>
std::function<R(Arg)> 
addFunctions(std::function<R(Arg)> first,
             std::function<R(Arg)> second) {
    return [=] (Arg a) {
        return first(a) + second(a);
    };
}

template<class R, class Arg0, class Arg1>
std::function<R(Arg0, Arg1)> 
addFunctions(std::function<R(Arg0, Arg1)> first,
             std::function<R(Arg0, Arg1)> second) {
    return [=] (Arg0 a0, Arg1 a1) {
        return first(a0, a1) + second(a0, a1);
    };
}

template<class R, class Arg0, class Arg1, class Arg2>
std::function<R(Arg0, Arg1, Arg2)> 
addFunctions(std::function<R(Arg0, Arg1, Arg2)> first,
             std::function<R(Arg0, Arg1, Arg2)> second) {
    return [=] (Arg0 a0, Arg1 a1, Arg2 a2) {
        return first(a0, a1, a2) + second(a0, a1, a2);
    };
}

}