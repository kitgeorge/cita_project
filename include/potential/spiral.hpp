#include <functional>
#include <cmath>
#include "force.hpp"
#include "potential_funcs.hpp"

namespace potential {

PotentialFuncs getSpiralPotential(double m, double k_R, 
                                  double amplitude,
                                  double pattern_speed);

}