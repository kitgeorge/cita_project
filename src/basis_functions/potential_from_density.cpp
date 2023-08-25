#include "potential_from_density.hpp"
#include <iostream>

namespace basis_functions {

template <typename DensityType>
PotentialFromDensity::
PotentialFromDensity(const DensityType& density, const BFE& expansion,
                     int n_max, int l_max): nl_max({n_max, l_max}),
                     coefficients(getCoefficients(density, expansion)),
                     trunc_density(getTruncDensity(expansion)),
                     trunc_potential(getTruncPotential(expansion)),
                     trunc_force(getTruncForce(expansion)) {}

PotentialFromDensity::
PotentialFromDensity(const std::vector<std::vector<std::complex<double>>> 
                     coefficients_, const BFE& expansion,
                     int n_max, int l_max): nl_max({n_max, l_max}),
                     coefficients(coefficients_),
                     trunc_density(getTruncDensity(expansion)),
                     trunc_potential(getTruncPotential(expansion)),
                     trunc_force(getTruncForce(expansion)) {}

PotentialFromDensity::
PotentialFromDensity(const PotentialFromDensity& old):
                     nl_max(old.nl_max), coefficients(old.coefficients),
                     trunc_density(old.trunc_density),
                     trunc_potential(old.trunc_potential),
                     trunc_force(old.trunc_force) {}

template <typename DensityType>
std::vector<std::vector<std::complex<double>>> 
PotentialFromDensity::
getCoefficients(const DensityType& density, const BFE& expansion) const {
    std::vector<std::vector<std::complex<double>>>
    output(nl_max[0] + 1, std::vector<std::complex<double>>(nl_max[1] + 1));
    std::vector<std::function<std::complex<double>()>>
    coefficient_functions((nl_max[0] + 1)*(nl_max[1] + 1));
    for(int i = 0; i <= nl_max[0]; ++i) {
        for(int j = 0; j <= nl_max[1]; ++j) {
            coefficient_functions[i*(nl_max[1] + 1) + j] = [=] () {
                std::cout << "Calculating BFE coefficients: "
                          << i << ", " << j << std::endl;
                return expansion.getCoefficient(i, j, density);
            };
        }
    }
    std::vector<std::complex<double>>
    output_flat = multithreading::executeInParallel(coefficient_functions);
    for(int i = 0; i <= nl_max[0]; ++i) {
        for(int j = 0; j <= nl_max[1]; ++j) {
            output[i][j] = output_flat[i*(nl_max[1] + 1) + j];
        }
    }
    return output;
}

template <typename DataType>
std::vector<std::vector<std::function<DataType(double, double)>>>
PotentialFromDensity::
getTerms(std::function<std::function<DataType(double, double)>
                         (int, int)> BFE_member_function,
                      const BFE& expansion) const {
    // Initialising table for 0 <= n <= n_max, -l_max <= 0 <= l_max
    std::vector<std::vector<std::function<DataType(double, double)>>>
    output(nl_max[0] + 1, 
            std::vector<std::function<DataType(double, double)>>(2*nl_max[1] + 1));
    // Populates table with respective basis functions multiplied by coefficients,
    // and complex conjugate for negative l (not sure if it should be a pointer
    // to expansion)
    for(int i = 0; i <= nl_max[0]; ++i) {
        for(int j = 0; j <= nl_max[1]; ++j) {
            output[i][nl_max[1] + j] = utility::multiplyFunction(BFE_member_function(i, j),
                                            coefficients[i][j]);
            if(j != 0) {
                output[i][nl_max[1] - j] 
                        = utility::conjugateFunction(output[i][nl_max[1] + j]);
            }
        }
    }
    return output;
}

// Trying to generalise the summation of basis functions for density, potential and
// force (in some ways this is quite ugly. DataType refers to complex or 
// complex array (force))
template <typename DataType>
std::function<DataType(double, double)>
PotentialFromDensity::
getTruncFunction(std::function<std::function<DataType(double, double)>
                    (int, int)> BFE_member_function,
                 const BFE& expansion) const {

    return utility::addFunctions(
                utility::flatten(getTerms(BFE_member_function,
                                          expansion)));
}

std::function<double(double, double)>
PotentialFromDensity::
getTruncDensity(const BFE& expansion) const {
    std::function<std::function<std::complex<double>(double, double)>(int, int)>
    wrapper = [=] (int i, int j) {
        return expansion.rho(i, j);
    };
    return utility::realFunction(getTruncFunction(wrapper, expansion));
}

std::function<double(double, double)>
PotentialFromDensity::
getTruncPotential(const BFE& expansion) const {
    std::function<std::function<std::complex<double>(double, double)>(int, int)>
    wrapper = [=] (int i, int j) {
        return expansion.psi(i, j);
    };
    return utility::realFunction(getTruncFunction(wrapper, expansion));
}

std::function<std::array<double, 2>(double, double)>
PotentialFromDensity::
getTruncForce(const BFE& expansion) const {
    std::function<std::function<std::array<std::complex<double>, 2>(double, double)>(int, int)>
    wrapper = [=] (int i, int j) {
        return expansion.psi_f(i, j);
    };
    return utility::realFunction(getTruncFunction(wrapper, expansion));
}

std::vector<std::vector<std::complex<double>>>
PotentialFromDensity::getCoefficients() const {
    return coefficients;
}

double PotentialFromDensity::calculateAbsNorm() const {
    double output = 0;
    for(int i = 0; i <= nl_max[0]; ++i) {
        for(int j = 0; j <= nl_max[1]; ++j) {
            output += pow(std::abs(coefficients[i][j]), 2);
        }
    }
    return sqrt(output);
}

template PotentialFromDensity::
         PotentialFromDensity(const std::function<double(double, double)>& density,
                              const BFE& expansion, int N_n, int N_l);
template PotentialFromDensity::
         PotentialFromDensity(const std::vector<std::array<double, 3>>& density,
                              const BFE& expansion, int N_n, int N_l);
template std::vector<std::vector<std::complex<double>>>
         PotentialFromDensity::
         getCoefficients(const std::function<double(double, double)>& density,
                         const BFE& expansion) const;
template std::vector<std::vector<std::complex<double>>>
         PotentialFromDensity::
         getCoefficients(const std::vector<std::array<double, 3>>& density,
                         const BFE& expansion) const;
template std::vector<std::vector<std::function<std::complex<double>(double, double)>>>
         PotentialFromDensity::
         getTerms(std::function<std::function<std::complex<double>(double, double)>
                    (int, int)> BFE_member_function,
                 const BFE& expansion) const;
template std::vector<std::vector<std::function<std::array<std::complex<double>, 2>(double, double)>>>
         PotentialFromDensity::
         getTerms(std::function<std::function<std::array<std::complex<double>, 2>(double, double)>
                    (int, int)> BFE_member_function,
                 const BFE& expansion) const;
template std::function<std::complex<double>(double, double)>
        PotentialFromDensity::
        getTruncFunction(std::function<std::function<std::complex<double>(double, double)>
                    (int, int)> BFE_member_function,
                 const BFE& expansion) const;
template std::function<std::array<std::complex<double>, 2>(double, double)>
        PotentialFromDensity::
        getTruncFunction(std::function<std::function<std::array<std::complex<double>, 2>(double, double)>
                    (int, int)> BFE_member_function,
                 const BFE& expansion) const;


}