#include "gtest/gtest.h"
#include "potential_funcs.hpp"
#include "axsym_funcs.hpp"
#include "mestel.hpp"
#include <cmath>
#include <numbers>

using namespace potential;

TEST(PotentialFuncsTest, CartesianForceWorks) {
    std::function<double(double, double, double)>
    potential = [] (double R, double phi, double t) {
        return -1/R + cos(phi);
    };
    std::function<std::array<double, 2>(double, double, double)> 
    force = [] (double R, double phi, double  t) {
        std::array<double, 2>
        output = {-1/pow(R, 2), sin(phi)};
        return output;
    };
    PotentialFuncs pot(potential, force);

    double R_max = 10;
    int N_R_values = 100;
    int N_phi_values = 360;
    
    for(int i = 0; i < N_R_values; ++i) {
        for(int j = 0; j < N_phi_values; ++j) {
            double R = (double)(i + 1)/N_R_values*R_max;
            double phi = (double)j/N_phi_values*2*std::numbers::pi;
            std::array<double, 2> cart = pot.cartesian_force(R*cos(phi), R*sin(phi), 0);
            std::array<double, 2> 
            cart_conv = {cart[0]*cos(phi) + cart[1]*sin(phi),
                         cart[1]*cos(phi) - cart[0]*sin(phi)};
            for(int i = 1; i < 2; ++i) {
                EXPECT_EQ(pow(pot.polar_force(R, phi, 0)[i] - cart_conv[i], 2) < 1e-6, 1);
            }
        }
    }



    // double R = 1;
    // double phi = std::numbers::pi/2;
    // double x = R*cos(phi);
    // double y = R*sin(phi);
    // std::array<double, 2> cart = pot.cartesian_force(x, y, 0);
    // std::array<double, 2> polar = {cart[0]*cos(phi) + cart[1]*sin(phi),
    //                                cart[1]*cos(phi) - cart[0]*sin(phi)};
    // for(int i = 0; i < 2; ++i) {
    //     ASSERT_EQ(pow(pot.polar_force(R, phi, 0)[i] - polar[i], 2) < 1e-6, 1);
    // }
}

TEST(PotentialFuncsTest, AddAxsymWorks) {
    AxsymFuncs as = getMestel(220*1000);
    double norm = as.potential(8*3e19, 0, 0);

    std::function<double(double, double, double)>
    potential = [norm] (double R, double phi, double t) {
        return norm/100*cos(phi);
    };
    std::function<std::array<double, 2>(double, double, double)> 
    force = [norm] (double R, double phi, double  t) {
        std::array<double, 2>
        output = {0, norm/100*sin(phi)/R};
        return output;
    };
    PotentialFuncs pot(potential, force);
    PotentialFuncs sum(pot, as);
    EXPECT_EQ(sum.potential(8*3e19, 0, 0), as.potential(8*3e19, 0, 0) + pot.potential(8*3e19, 0, 0));
    for(int i = 0; i < 2; ++i) {
        EXPECT_EQ(sum.polar_force(8*3e19, 0, 0)[i], as.polar_force(8*3e19, 0, 0)[i] + pot.polar_force(8*3e19, 0, 0)[i]);
        EXPECT_EQ(sum.cartesian_force(8*3e19, 0, 0)[i], as.cartesian_force(8*3e19, 0, 0)[i] + pot.cartesian_force(8*3e19, 0, 0)[i]);
    }
}