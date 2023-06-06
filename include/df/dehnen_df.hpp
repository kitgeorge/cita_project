#include <cmath>
#include <functional>
#include "axsym_funcs.hpp"

namespace df {

std::function<double(double, double)> 
getDehnenDF(const potential::AxsymFuncs& pot, 
            std::function<double(double)> surface_density,
            std::function<double(double)> sigma_R);

}