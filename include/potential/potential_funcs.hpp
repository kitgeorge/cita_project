#pragma once
#include "force.hpp"
#include "coords.hpp"
#include "add_functions.hpp"
#include "potential_from_density.hpp"
#include <array>
#include <assert.h>

namespace vrs = vectors;

namespace potential {

struct PotentialFuncs {
    double potential(double R, double phi, double t) const;
    std::array<double, 2> polar_force(double R, double phi, double t) const;
    std::array<double, 2> cartesian_force(double x, double y, double t) const;




    std::vector<std::function<double(double, double, double)>> potentials;
    std::vector<std::function<std::array<double, 2>(double, double, double)>>
    polar_forces;
    std::vector<std::function<std::array<double, 2>(double, double, double)>>
    cartesian_forces;

    PotentialFuncs(std::function<double(double, double, double)> potential_,
                   std::function<std::array<double, 2>(double, double, double)> 
                   polar_force_);

    PotentialFuncs(const basis_functions::PotentialFromDensity& p);

    PotentialFuncs(const std::vector<PotentialFuncs>& p);

    PotentialFuncs(const PotentialFuncs& old);

    void operator +=(const PotentialFuncs& p);
    PotentialFuncs operator *(const double& factor) const;

    void multiply(std::function<double(double, double, double)> envelope);
    private:
        std::function<std::array<double, 2>(double, double, double)>
        getCartesianForce(std::function<std::array<double, 2>(double, double, double)>
                          polar_force);
};

}