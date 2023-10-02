#include "axsym_funcs.hpp"

namespace potential {

double AxsymFuncs::EGivenPolar(std::array<std::array<double, 2>, 2> coords) {
    double T = 0.5*(pow(coords[1][0], 2) + pow(coords[1][1], 2));
    double V = potential_R(coords[0][0]);
    return T + V;
}

double AxsymFuncs::LGivenPolar(std::array<std::array<double, 2>, 2> coords) {
    return coords[0][0]*coords[1][1];
}

}