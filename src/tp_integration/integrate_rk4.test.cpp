#include "gtest/gtest.h"
#include "integrate_rk4.hpp"
#include <iostream>

using namespace tp_integration;
namespace vrs = vectors;
namespace ptl = potential;

TEST(Rk4Test, HarmonicOscillator) {
    double kappa = 1;
    std::function<std::array<double, 2>(double, double, double)> 
    force = [kappa](double R, double phi, double t) {
        std::array<double, 2> output = {-kappa*R, 0};
        return output;
    };
    std::function<double(double, double, double)>
    potential = [kappa](double R, double phi, double t) {
        return kappa/2*pow(R, 2);
    };
    ptl::PotentialFuncs pot(potential, force);
    vrs::Coords2d initial_conditions({{{{1, 0}}, {{0, 0}}}}, 0);
    double timestep = 0.01;
    int N_timesteps = 1000;
    std::vector<vrs::Coords2d>
    trajectory = getTpIntegrationFunction(pot, initial_conditions,
                                          0, timestep, N_timesteps)();
    std::vector<vrs::Coords2d>
    trajectory2 = getTpIntegrationFunction(pot, initial_conditions,
                                          0, 10*timestep, N_timesteps)();
    std::function<double(vrs::Coords2d)>
    energy = [potential] (vrs::Coords2d coords) {
        double output = (pow(coords.cartesian[1][0], 2) 
                         + pow(coords.cartesian[1][1], 2))/2;
        output += potential(coords.polar[0][0], coords.polar[0][1], 0);
        return output;
    };
    // std::cout << energy(trajectory[N_timesteps - 1]) - energy(trajectory[0])
    //           << ", " 
    //           << energy(trajectory2[N_timesteps - 1]) - energy(trajectory2[0])
    //           << std::endl; 
    ASSERT_EQ(energy(trajectory[N_timesteps - 1]) - energy(trajectory[0]) 
                    < 1e-5*energy(trajectory[0]), 1);
}