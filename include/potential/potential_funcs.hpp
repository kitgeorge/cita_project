#pragma once
#include "force.hpp"
#include "coords.hpp"
#include "add_functions.hpp"
#include <array>

namespace vrs = vectors;

namespace potential {

struct PotentialFuncs {
    std::function<double(double, double, double)> potential;
    std::function<std::array<double, 2>(double, double, double)> polar_force;
    std::function<std::array<double, 2>(double, double, double)> cartesian_force;

    PotentialFuncs(std::function<double(double, double, double)> potential_,
                   std::function<std::array<double, 2>(double, double, double)> 
                   polar_force_):
                   potential(potential_), polar_force(polar_force_) {
        cartesian_force = [this](double x, double y, double t) {
            std::array<std::array<double, 2>, 2>
            coords = {{ {{x, y}}, {{0, 0}} }};
            vrs::Coords2d pos(coords, 0);
            std::array<double, 2> f = {polar_force(pos.polar[0][0], 
                                                   pos.polar[0][1], t)[0],
                                       polar_force(pos.polar[0][0], 
                                                   pos.polar[0][1], t)[1]};
            f = vrs::getCartesianVector2d(pos.polar[0], f);
            return f;
        };
    }
    PotentialFuncs(const PotentialFuncs& p0, const PotentialFuncs& p1) {
        potential = utility::addFunctions(p0.potential, p1.potential);
        polar_force = utility::addFunctions(p0.polar_force, p1.polar_force);
        cartesian_force = utility::addFunctions(p0.cartesian_force,
                                                p1.cartesian_force);
    }
    PotentialFuncs(const PotentialFuncs& old) {
        potential = old.potential;
        polar_force = old.polar_force;
        cartesian_force = old.cartesian_force;
    }
};

}