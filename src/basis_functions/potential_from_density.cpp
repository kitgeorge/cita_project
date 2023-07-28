#include "potential_from_density.hpp"

namespace basis_functions {

template <typename DensityType>
PotentialFromDensity::
PotentialFromDensity(const DensityType& density, const BFE& expansion,
                     int n_max, int l_max): nl_max({n_max, l_max}),
                     coefficients(getCoefficients(density, expansion)),
                     trunc_density(getTruncDensity(expansion)),
                     trunc_potential(getTruncPotential(expansion)),
                     trunc_force(getTruncForce(expansion)) {}

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

// Trying to generalise the summation of basis functions for density, potential and
// force (in some ways this is quite ugly. DataType refers to complex or 
// complex array (force))
template <typename DataType>
std::function<DataType(double, double)>
PotentialFromDensity::
getTruncFunction(std::function<std::function<DataType(double, double)>
                    (int, int)> BFE_member_function,
                 const BFE& expansion) const {

    // Initialising table for 0 <= n <= n_max, -l_max <= 0 <= l_max
    std::vector<std::vector<std::function<DataType(double, double)>>>
    terms(nl_max[0] + 1, 
          std::vector<std::function<DataType
                (double, double)>>(2*nl_max[1] + 1)); // Includes negative l
    
    // Populates table with respective basis functions multiplied by coefficients,
    // and complex conjugate for negative l (not sure if it should be a pointer
    // to expansion)
    for(int i = 0; i <= nl_max[0]; ++i) {
        for(int j = 0; j <= nl_max[1]; ++j) {
            terms[i][nl_max[1] + j] 
                = utility::multiplyFunction(BFE_member_function(i, j), 
                                            coefficients[i][j]);
            if(j != 0) {
                terms[i][nl_max[1] - j]
                    = utility::conjugateFunction(terms[i][nl_max[1] + j]);
            }
        }
    }
    
    return utility::addFunctions(utility::flatten(terms));
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