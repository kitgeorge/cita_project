#include "add_functions.hpp"

namespace utility {

double addFunctions(const std::vector<std::function<double(double, double, double)>>& 
                    functions, double a0, double a1, double a2) {
    double output = 0;
    for(auto f : functions) {
        output += f(a0, a1, a2);
    }
    return output;
}

}