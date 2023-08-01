#include "bfe.hpp"
#include "potential_from_density.hpp"
#include "vector_io.hpp"
#include "flatten.hpp"
#include "comp_funcs.hpp"
#include "mestel_spiral.hpp"
#include "execute_in_parallel.hpp"
#include "gamma.hpp"
#include "nd_vectors.hpp"
#include <iostream>


int main() {

    int n_max = 30;
    int l_max = 30;
    int k_Ka = 4;
    double R_Ka = 20*Units::kpc;
    int N_R = 1000;
    int N_phi = 1000;

    basis_functions::BFE expansion(k_Ka, R_Ka, N_R, N_phi);

    std::function<double(double, double, double)>
    spiral_density_time = potential::getStandardSpiralDensity();
    std::function<double(double, double)>
    spiral_density = [spiral_density_time] (double R, double phi) {
        return spiral_density_time(R, phi, 0);
    };
    basis_functions::PotentialFromDensity
    pot(spiral_density, expansion, n_max, l_max);

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

    for(int i = 0; i < N_R_values; ++i) {
        double R = (double)i/N_R_values*R_Ka;
        for(int j = 0; j < N_phi_values; ++j) {
            double phi = (double)j/N_phi_values*2*std::numbers::pi;
            density_functions[i*N_phi_values + j] = [=] () {
                return spiral_density(R, phi);
            };
            trunc_density_functions[i*N_phi_values + j] = [&pot, R, phi] () {
                return pot.trunc_density(R, phi);
            };
            trunc_potential_functions[i*N_phi_values + j] = [&pot, R, phi] () {
                return pot.trunc_potential(R, phi);
            };
            trunc_force_functions[i*N_phi_values + j] = [&pot, R, phi] () {
                return pot.trunc_force(R, phi);
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
    utility::writeCsv("../data/basis_functions/spiral_trunc_density.csv",
                      trunc_density_values);
    utility::writeCsv("../data/basis_functions/spiral_trunc_potential.csv",
                      trunc_potential_values);
    utility::writeCsv("../data/basis_functions/spiral_trunc_force.csv",
                      utility::flatten(trunc_force_values));
    utility::writeCsv("../data/basis_functions/spiral_density.csv",
                      density_values);

    std::array<int, 3> coefficients_shape = {n_max + 1, l_max + 1, 2};
    utility::vector3d<double> 
    coefficients = utility::makeShape<double>(coefficients_shape);
    for(int i = 0; i <= n_max; ++i) {
        for(int j = 0; j <= l_max; ++j) {
            coefficients[i][j][0] = std::abs(pot.getCoefficients()[i][j]);
            coefficients[i][j][1] = std::arg(pot.getCoefficients()[i][j]);
        }
    }
    utility::writeCsv("../data/basis_functions/spiral_bfe_coefficients.csv",
                      utility::flatten(coefficients));
}
