#include <array>
#include <cstdlib>
#include <functional>
#include <numbers>

namespace df {

std::array<std::array<double, 2>, 2> 
getDFSample(std::function<double(double, double)> df, 
                          double R_min, double R_max, 
                          double v_R_max, double delta_v_phi_max, double v_c);


std::array<double, 2>
getDFSampleEL(std::function<double(double, double)> df,
              std::array<std::array<double, 2>, 2> bounds);
}