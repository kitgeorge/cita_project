#include <functional>
#include <numbers>
#include <cmath>
#include <array>
#include "axsym_funcs.hpp"
#include "mapXVtoAA2D.hpp"
#include "units.hpp"

namespace actions {

std::function<double(double)> 
getThetaRSTQIntegrand(potential::AxsymFuncs pot,
                      double E, double L);

std::array<double, 2>
findApsis(potential::AxsymFuncs pot, double E, double L);

std::function<double(double)>
getPhiEff(potential::AxsymFuncs pot, double L);

}