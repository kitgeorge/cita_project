#include "bfe.hpp"
#include "potential_from_density.hpp"
#include "gtest/gtest.h"
#include "vector_io.hpp"
#include "flatten.hpp"
#include "comp_funcs.hpp"
#include "mestel_spiral.hpp"
#include "execute_in_parallel.hpp"
#include "comp_funcs.hpp"
#include <optional>

using namespace basis_functions;

class BFETest : public ::testing::Test {
    protected:
        std::optional<BFE> expansion;
        int n_max;
        int l_max;
        double R_Ka;
        int N_R;
        int N_phi;

        double mestel_vc;
        double mestel_R0;
        
        std::optional<PotentialFromDensity> pot;

        std::function<double(double, double)> mestel_density;
        std::function<double(double, double)> mestel_potential;
        std::function<std::array<double, 2>(double, double)>
        mestel_force;

        void SetUp() override {
            R_Ka = 20*Units::kpc;
            N_R = 5000;
            N_phi = 5000;
            expansion.emplace(R_Ka, N_R, N_phi);
            // N_R and N_phi are being used for too many things
            // but can't be bothered to fix yet
            n_max = 30;
            l_max = 30;
            mestel_vc = 220*Units::kms;
            mestel_R0 = 8*Units::kpc;
            double vc = mestel_vc;
            double R0 = mestel_R0;
            mestel_density = [vc] (double R, double phi) {
                return pow(vc, 2)/(2*std::numbers::pi*Units::G*R);
            };
            // pot.emplace(mestel_density, 
            //             expansion.value(), n_max, l_max);
            mestel_potential = [vc, R0] (double R, double phi) {
                return pow(vc, 2)*log(R/R0);
            };
            mestel_force = [vc] (double R, double phi) {
                std::array<double, 2> 
                output = {-pow(vc, 2)/R, 0};
                return output;
            };
        }
};

TEST_F(BFETest, DISABLED_CheckOrthonormal) {
    // Very coarse check as unit test. Integrals with more intervals
    // are expected to converge to orthonormality, but this takes longer 
    // to check so we won't do that every time
    for(int h = 0; h <= n_max; ++h) {
        for(int i = 0; i < n_max; ++i) {
            for(int j = 0; j <= l_max; ++j) {
                for(int k = 0; k <= l_max; ++k) {
                    // std::cout << h << ", " << i << ", " << j << ", " << k << std::endl;
                    std::function<std::complex<double>(double, double)>
                    conj = [this, h, j] (double R, double theta) {
                        return std::conj(expansion.value().psi(h, j)(R, theta));
                    };
                    if(h == i && j == k) {
                        EXPECT_LT(std::abs(scalarProduct(
                                                conj, expansion.value().rho(i, k),
                                                R_Ka, N_R, N_phi) + 1.0),
                                                0.01);
                    }
                    else {
                        EXPECT_LT(std::abs(scalarProduct(
                                                conj, expansion.value().rho(i, k),
                                                R_Ka, N_R, N_phi)),
                                                0.01);
                    }
                }
            }
        }
    }
}

TEST_F(BFETest, DISABLED_CheckNorm) {
    for(int i = 0; i <= n_max; ++i) {
        std::cout << i << std::endl;
        EXPECT_LT(std::abs(scalarProduct(utility::conjugateFunction(expansion.value().psi(i, 2)),
                                         expansion.value().rho(i, 2), R_Ka, N_R, N_phi) + 1.0), 0.25);
    }
}

TEST_F(BFETest, DISABLED_CheckCoefficient) {
    for(int h = 0; h < n_max; ++h) {
        for(int i = 0; i < n_max; ++i) {
            for(int j = 0; j < l_max; ++j) {
                for(int k = 0; k < l_max; ++k) {
                    // std::cout << h << ", " << i << ", " << j << ", " << k << std::endl;
                    std::vector<std::function<std::complex<double>(double, double)>>
                    realBasisFunc = {expansion.value().rho(i, k),
                                             utility::conjugateFunction(expansion.value().rho(i, k))};
                    std::complex<double>
                    coefficient = expansion.value().getCoefficient(h, j, 
                                        utility::realFunction((utility::addFunctions(realBasisFunc))));
                    if(k == 0) {
                        coefficient = coefficient/2.0;
                    }
                    if(h == i && j == k) {
                        EXPECT_LT(std::abs(coefficient.real() - 1), 0.25);
                    }
                    else {
                        EXPECT_LT(std::abs(coefficient), 0.25);
                    }
                }
            }
        }
    }
}

TEST_F(BFETest, DISABLED_PlotBasisFunctions) {
    int N_n = 6;
    int l = 2;
    int N_k = 6;
    std::vector<double> 
    basis_function_values(N_n*2*N_R);
    for(int h = 0; h < N_n; ++h) {
        BFE funcs(R_Ka, N_R, N_phi);
        std::array<std::function<std::complex<double>(double, double)>, 2>
        basis_functions = {funcs.psi(h, l), funcs.rho(h, l)};
        for(int j = 0; j < 2; ++j) {
            for(int k = 0; k < N_R; ++k) {
                double R = (double)k/N_R*R_Ka;
                basis_function_values[h*2*N_R + j*N_R + k]
                    = basis_functions[j](R, 0).real();
            }
        }
    }
    utility::writeCsv("../test_data/basis_functions/basis_functions.csv",
                      basis_function_values);
}

TEST_F(BFETest, DISABLED_CheckMestel) {
    // Plot truncated functions against R (phi = 0)
    int N_R_values = 1000;


    std::vector<std::function<double()>> trunc_density_functions(N_R_values);
    std::vector<std::function<double()>> trunc_potential_functions(N_R_values);
    std::vector<std::function<std::array<double, 2>()>> 
    trunc_force_functions(N_R_values);

    std::vector<std::function<double()>> density_functions(N_R_values);
    std::vector<std::function<double()>> potential_functions(N_R_values);
    std::vector<std::function<std::array<double, 2>()>> 
    force_functions(N_R_values);

    for(int i = 0; i < N_R_values; ++i) {
        double R = (double)i/N_R_values*R_Ka;
        trunc_density_functions[i] = [R, this] () {
            return pot.value().trunc_density(R, 0);
        };
        trunc_potential_functions[i] = [R, this] () {
            return pot.value().trunc_potential(R, 0);
        };
        trunc_force_functions[i] = [R, this] () {
            return pot.value().trunc_force(R, 0);
        };
        density_functions[i] = [R, this] () {
            return mestel_density(R, 0);
        };
        potential_functions[i] = [R, this] () {
            return mestel_potential(R, 0);
        };
        force_functions[i] = [R, this] () {
            return mestel_force(R, 0);
        };
    };
    std::cout << "A" << std::endl;
    std::vector<double> 
    trunc_density_values = multithreading::
                           executeInParallel(trunc_density_functions);
    std::cout << "A" << std::endl;
    std::vector<double> 
    trunc_potential_values = multithreading::
                             executeInParallel(trunc_potential_functions);
    std::cout << "A" << std::endl;
    std::vector<std::array<double, 2>> 
    trunc_force_values = multithreading::
                         executeInParallel(trunc_force_functions);

    std::cout << "A" << std::endl;
    std::vector<double> 
    density_values = multithreading::
                           executeInParallel(density_functions);
    std::cout << "A" << std::endl;
    std::vector<double> 
    potential_values = multithreading::
                             executeInParallel(potential_functions);
    std::cout << "A" << std::endl;
    std::vector<std::array<double, 2>> 
    force_values = multithreading::
                         executeInParallel(force_functions);


    utility::writeCsv("../test_data/basis_functions/trunc_density.csv",
                      trunc_density_values);
    utility::writeCsv("../test_data/basis_functions/trunc_potential.csv",
                      trunc_potential_values);
    utility::writeCsv("../test_data/basis_functions/trunc_force.csv",
                      utility::flatten(trunc_force_values));

    utility::writeCsv("../test_data/basis_functions/mestel_density.csv",
                      density_values);
    utility::writeCsv("../test_data/basis_functions/mestel_potential.csv",
                      potential_values);
    utility::writeCsv("../test_data/basis_functions/mestel_force.csv",
                      utility::flatten(force_values));
}

TEST_F(BFETest, DISABLED_CheckSpiral) {
    std::function<double(double, double, double)>
    spiral_density_time = potential::getStandardSpiralDensity();
    std::function<double(double, double)>
    spiral_density = [spiral_density_time] (double R, double phi) {
        return spiral_density_time(R, phi, 0);
    };
    std::cout << spiral_density(8, 0) << std::endl;
    int N_R_values = 100;
    int N_phi_values = 360;
    std::vector<std::function<double()>>
    density_functions(N_R_values*N_phi_values);
    std::vector<std::function<double()>> 
    trunc_density_functions(N_R_values*N_phi_values);
    std::vector<std::function<double()>> 
    trunc_potential_functions(N_R_values*N_phi_values);
    std::vector<std::function<std::array<double, 2>()>> 
    trunc_force_functions(N_R_values*N_phi_values);

    PotentialFromDensity p(std::make_shared<BFE>(expansion.value()), n_max, l_max);
    p.initFromDensity(spiral_density);

    for(int i = 0; i < N_R_values; ++i) {
        double R = (double)i/N_R_values*R_Ka;
        for(int j = 0; j < N_phi_values; ++j) {
            double phi = (double)j/N_phi_values*2*std::numbers::pi;
            density_functions[i*N_phi_values + j] = [=] () {
                return spiral_density(R, phi);
            };
            trunc_density_functions[i*N_phi_values + j] = [&p, R, phi] () {
                return p.trunc_density(R, phi);
            };
            trunc_potential_functions[i*N_phi_values + j] = [&p, R, phi] () {
                return p.trunc_potential(R, phi);
            };
            trunc_force_functions[i*N_phi_values + j] = [&p, R, phi] () {
                return p.trunc_force(R, phi);
            };
        }
    }
    std::vector<double> 
    density_values = multithreading::executeInParallel(density_functions);
    std::vector<double>
    trunc_density_values = multithreading::
                                executeInParallel(trunc_density_functions);
    std::vector<double>
    trunc_potential_values = multithreading::
                        executeInParallel(trunc_potential_functions);
    std::vector<std::array<double, 2>>
    trunc_force_values = multithreading::
                        executeInParallel(trunc_force_functions);
    utility::writeCsv("../test_data/basis_functions/spiral_trunc_density.csv",
                      trunc_density_values);
    utility::writeCsv("../test_data/basis_functions/spiral_trunc_potential.csv",
                      trunc_potential_values);
    utility::writeCsv("../test_data/basis_functions/spiral_trunc_force.csv",
                      utility::flatten(trunc_force_values));
    utility::writeCsv("../test_data/basis_functions/spiral_density.csv",
                      density_values);
}