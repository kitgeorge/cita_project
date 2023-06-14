#include "gtest/gtest.h"
#include "spiral.hpp"
#include <numbers>
#include "vector_io.hpp"
#include <cmath>

using potential::getSpiralPotential;
namespace vrs = vectors;

TEST(GetSpiralPotentialTest, CartForceConsistent) {
    potential::PotentialFuncs 
    p = getSpiralPotential(2, 1, 1, 2*std::numbers::pi, 0);
    double R_max = 10;
    double t_max = 1;
    int N_R_values = 100;
    int N_phi_values = 360;
    int N_t_values = 10;

    for(int i = 0; i < N_R_values; ++i) {
        for(int j = 0; j < N_phi_values; ++j) {
            for(int k = 0; k < N_t_values; ++k) {
                double R = (double)(i + 1)/N_R_values*R_max;
                double phi = (double)j/N_phi_values*2*std::numbers::pi;
                double t = (double)k/N_t_values*t_max;
                std::array<double, 2> cart = p.cartesian_force(R*cos(phi), R*sin(phi), t);
                std::array<double, 2> 
                cart_conv = {cart[0]*cos(phi) + cart[1]*sin(phi),
                             cart[1]*cos(phi) - cart[0]*sin(phi)};
                for(int i = 0; i < 2; ++i) {
                    EXPECT_EQ(pow(p.polar_force(R, phi, t)[i] - cart_conv[i], 2) < 1e-6, 1);
                }
            }
        }
    }
}

TEST(GetSpiralPotentialTest, PlottingCheck) {
    potential::PotentialFuncs 
    p = getSpiralPotential(2, 1, 1, 2*std::numbers::pi, 0);
    double R_max = 10;
    double t_max = 1;
    int N_R_values = 100;
    int N_phi_values = 360;
    int N_t_values = 10;

    std::vector<double> potential_data(N_R_values*N_phi_values*N_t_values);
    std::vector<double>
    force_data(N_R_values*N_phi_values*N_t_values*2);
    std::vector<double>
    cart_force_data(N_R_values*N_phi_values*N_t_values*2); 


    for(int i = 0; i < N_R_values; ++i) {
        for(int j = 0; j < N_phi_values; ++j) {
            for(int k = 0; k < N_t_values; ++k) {
                double R = (double)(i + 1)/N_R_values*R_max;
                double phi = (double)j/N_phi_values*2*std::numbers::pi;
                double t = (double)k/N_t_values*t_max;
                potential_data[i*N_phi_values*N_t_values
                               + j*N_t_values + k] = p.potential(R, phi, t);
                force_data[i*N_phi_values*N_t_values*2
                           + j*N_t_values*2 + k*2] = p.polar_force(R, phi, t)[0];
                force_data[i*N_phi_values*N_t_values*2
                           + j*N_t_values*2 + k*2 + 1] 
                    = p.polar_force(R, phi, t)[1];
                cart_force_data[i*N_phi_values*N_t_values*2
                           + j*N_t_values*2 + k*2] 
                    = p.cartesian_force(R*cos(phi), R*sin(phi), t)[0];
                cart_force_data[i*N_phi_values*N_t_values*2
                           + j*N_t_values*2 + k*2 + 1] 
                    = p.cartesian_force(R*cos(phi), R*sin(phi), t)[1];

                // std::array<double, 2> cart = {cart_force_data[i*N_phi_values*N_t_values*2
                //            + j*N_t_values*2 + k*2], cart_force_data[i*N_phi_values*N_t_values*2
                //            + j*N_t_values*2 + k*2 + 1]  };
                // std::array<double, 2>
                // cart_conv = {cart[0]*cos(phi) + cart[1]*sin(phi),
                //              cart[1]*cos(phi) - cart[0]*sin(phi)};
                // for(int l = 0; l < 2; ++l) {
                //     EXPECT_EQ(pow(force_data[i*N_phi_values*N_t_values*2 
                //                          + j*N_t_values*2 + k*2 + l] -
                //               cart_conv[l], 2)< 1e-6, 1);
                // }
            }
        }
    }
    utility::writeCsv("../test_data/potential/spiral_potential.csv",
                     potential_data);
    utility::writeCsv("../test_data/potential/spiral_force.csv",
                     force_data);
    utility::writeCsv("../test_data/potential/spiral_force_cart.csv",
                      cart_force_data);
}