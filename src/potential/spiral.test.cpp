#include "gtest/gtest.h"
#include "spiral.hpp"
#include <numbers>
#include "vector_io.hpp"

using potential::getSpiralPotential;
namespace vrs = vectors;

TEST(GetSpiralPotentialTest, PlottingCheck) {
    potential::PotentialFuncs 
    p = getSpiralPotential(2, 1, 1, 2*std::numbers::pi);
    double R_max = 10;
    double t_max = 1;
    int N_R_values = 100;
    int N_phi_values = 360;
    int N_t_values = 10;

    std::vector<double> potential_data(N_R_values*N_phi_values*N_t_values);
    std::vector<double>
    force_data(N_R_values*N_phi_values*N_t_values*2);


    for(int i = 0; i < N_R_values; ++i) {
        for(int j = 0; j < N_phi_values; ++j) {
            for(int k = 0; k < N_t_values; ++k) {
                double R = (double)(i + 1)/N_R_values*R_max;
                double phi = (double)j/N_phi_values*2*std::numbers::pi;
                double t = (double)k/N_t_values*t_max;
                potential_data[i*N_phi_values*N_t_values
                               + j*N_t_values + k] = p.potential(R, phi, t);
                force_data[i*N_phi_values*N_t_values*2
                           + j*N_t_values*2 + k*2] = p.force(R, phi, t).f_R;
                force_data[i*N_phi_values*N_t_values*2
                           + j*N_t_values*2 + k*2 + 1] 
                    = p.force(R, phi, t).f_phi;
            }
        }
    }
    utility::writeCsv("../test_data/potential/spiral_potential.csv",
                     potential_data);
    utility::writeCsv("../test_data/potential/spiral_force.csv",
                     force_data);
}