#include "force.hpp"
#include "add_functions.hpp"
namespace vrs = vectors;

namespace potential {

struct PotentialFuncs {
    std::function<double(double, double, double)> potential;
    std::function<vrs::Force(double, double, double)> force;

    PotentialFuncs(std::function<double(double, double, double)> potential_,
                   std::function<vrs::Force(double, double, double)> force_):
                   potential(potential_), force(force_) {};
    PotentialFuncs(PotentialFuncs p0, PotentialFuncs p1) {
        potential = utility::addFunctions(p0.potential, p1.potential);
        force = utility::addFunctions(p0.force, p1.force);
    }
};

}