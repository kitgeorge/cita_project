#include "gtest/gtest.h"
#include "potential_funcs.hpp"
#include "axsym_funcs.hpp"
#include "mestel.hpp"
#include "units.hpp"
#include "spiral.hpp"
#include "shape.hpp"
#include <cmath>
#include <numbers>
#include <iostream>

using namespace potential;

TEST(PotentialFuncsTest, DISABLED_PlottingPotentialFromDensity) {
    std::array<int, 2> shape = {65, 33};
    utility::vector2d<std::complex<double>>
    coefficients = utility::makeShape<std::complex<double>>(shape);
    for(int i = 0; i <= 64; ++i) {
        for(int j = 0; j <= 32; ++j) {
            coefficients[i][j] = 0;
        }
    }
    coefficients[0][1] = 1;
    basis_functions::PotentialFromDensity bfe_pot;
    bfe_pot.initFromCoefficients(coefficients);
    PotentialFuncs pot(bfe_pot);
    int N_R = 100;
    int N_phi = 300;
    std::array<int, 2> plot_shape = {N_R, N_phi};
    double R_max = 1;
    utility::vector2d<double> 
    potential_values = utility::makeShape<double>(plot_shape);
    utility::vector2d<std::array<double, 2>> 
    force_values = utility::makeShape<std::array<double, 2>>(plot_shape);
    for(int i = 0; i < plot_shape[0]; ++i) {
        for(int j = 0; j < plot_shape[1]; ++j) {
            double R = (double)i/plot_shape[0]*R_max;
            double phi = (double)j/plot_shape[1]*2*std::numbers::pi;
            potential_values[i][j] = pot.potential(R, phi, 0);
            force_values[i][j] = pot.polar_force(R, phi, 0);
        }
    }
    utility::writeCsv("../test_data/potential/pfd_potential.csv",
                      utility::flatten(potential_values));
    utility::writeCsv("../test_data/potential/pfd_polar_force.csv",
                      utility::flatten(utility::flatten(force_values)));

}

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
    AxsymFuncs as = getMestel(220*1000, 8);
    double norm = as.potential(8, 0, 0);

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
    PotentialFuncs sum({pot, as});
    EXPECT_EQ(sum.potential(8*3e19, 0, 0), as.potential(8*3e19, 0, 0) + pot.potential(8*3e19, 0, 0));
    for(int i = 0; i < 2; ++i) {
        EXPECT_EQ(sum.polar_force(8*3e19, 0, 0)[i], as.polar_force(8*3e19, 0, 0)[i] + pot.polar_force(8*3e19, 0, 0)[i]);
        EXPECT_EQ(sum.cartesian_force(8*3e19, 0, 0)[i], as.cartesian_force(8*3e19, 0, 0)[i] + pot.cartesian_force(8*3e19, 0, 0)[i]);
    }
}

TEST(PotentialFuncsTest, MultiplySpiralWorks) {
    double R_0 = 8*Units::kpc;
    AxsymFuncs as = getMestel(220*Units::kms, R_0);
    double spiral_amplitude_fraction = 0.001;
    double local_Omega = as.Omega(R_0);
    int m = 2;
    double k_R = 2*std::numbers::pi/Units::kpc;
    PotentialFuncs spiral = getSpiralPotential(2, k_R, spiral_amplitude_fraction/k_R
                                                       *(-1)*as.polar_force(R_0, 0, 0)[0], local_Omega, 0);
    std::function<double(double, double, double)>
    envelope = [R_0] (double R, double phi, double t) {
        return exp(-pow(R - R_0, 2)/(2*Units::kpc))*pow(cos(t/Units::Myr), 2);
    };
    PotentialFuncs spiral2 = spiral;
    spiral2.multiply(envelope);
    double R = R_0 + 0.7*Units::kpc;
    double phi = 3.2;
    double t = 0.56*Units::Myr;
    EXPECT_EQ(spiral2.potential(R, phi, t), envelope(R, phi, t)*spiral.potential(R, phi, t));
    for(int i = 0; i < 2; ++i) {
        EXPECT_EQ(spiral2.polar_force(R, phi, t)[i], envelope(R, phi, t)*spiral.polar_force(R, phi, t)[i]);
        EXPECT_EQ(spiral2.cartesian_force(R*cos(phi), R*sin(phi), t)[i], envelope(R, phi, t)*spiral.cartesian_force(R*cos(phi), R*sin(phi), t)[i]);
    }

}

TEST(PotentialFuncsTest, MestelResonanceWorks) {
    double v_c = 220*Units::kms;
    double R_0 = 8*Units::kpc;
    AxsymFuncs pot = getMestel(v_c, R_0);
    double R = pot.resonanceRadius(1, 2, pot.Omega(R_0));
    std::cout << R << std::endl;
}