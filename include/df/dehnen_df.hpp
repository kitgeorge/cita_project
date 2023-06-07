#include <cmath>
#include <functional>
#include <algorithm>
#include "axsym_funcs.hpp"

namespace df {

std::function<double(double, double)> 
getDehnenDF(const potential::AxsymFuncs& pot, 
            std::function<double(double)> surface_density,
            std::function<double(double)> sigma_R);

std::function<double(double, double)>
getNormDehnenDF(const potential::AxsymFuncs& pot,
                std::function<double(double)> surface_density,
                std::function<double(double)> sigma_R,
                std::array<double, 2> E_bounds);

}