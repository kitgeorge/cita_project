#pragma once
#include "bfe.hpp"
#include "add_functions.hpp"
#include "comp_funcs.hpp"
#include "flatten.hpp"
#include "execute_in_parallel.hpp"

namespace basis_functions {

class PotentialFromDensity {
    const basis_functions::BFE expansion;
    const std::array<int, 2> nl_max;
    std::vector<std::vector<std::complex<double>>> coefficients;


    template <typename DensityType>
    std::vector<std::vector<std::complex<double>>> 
    getCoefficients(const DensityType& density) const;

    // BFE_member_function refers to which set of basis functions
    // (psi, psi_f or rho) to expand in
    template <typename DataType>
    std::vector<std::vector<std::function<DataType(double, double)>>>
    getTerms(std::function<std::function<DataType(double, double)>
                (int, int)> BFE_member_function) const;

    // Sums the BFE for potential, force or density. I keep the above
    // function separate for diagnostic purposes
    template <typename DataType>
    std::function<DataType(double, double)>
    getTruncFunction(std::function<std::function<DataType(double, double)>
                        (int, int)> BFE_member_function) const;

    // Applies the above for each 
    std::function<double(double, double)> 
    getTruncDensity() const;
    std::function<double(double, double)> 
    getTruncPotential() const;
    std::function<std::array<double, 2>(double, double)> 
    getTruncForce() const;

    public:
        std::function<double(double, double)> trunc_density;
        std::function<double(double, double)> trunc_potential;
        std::function<std::array<double, 2>(double, double)> trunc_force;

        PotentialFromDensity(int n_max=64, int l_max=32);
        PotentialFromDensity(const basis_functions::BFE& expansion_,
                             int n_max=64, int l_max=32);
        PotentialFromDensity(const PotentialFromDensity& old);


        /**
         * Finds the mass of particles for a df to match a target density.
         * 
         * We start with a DF such as from Dehnen (1999), which approximates
         * a target density profile such as a Mestel disc. It is not trivial
         * to set the masses of particles (for the given number sampled) which
         * best approximate this density profile. We achieve that here by
         * approximating the smooth density of the DF with the BFE. We 
         * choose that every particle has the same mass, we set the mass to 1
         * and we compare the resulting truncated density
         * to a target density profile, at a given position. The correct mass
         * is the factor by which target density is greater.
         *
         * @param positions Vector of positions of stars, {{R, phi}}
         * @param target_density target density function, rho(R, phi)
         * @param target_coords {R, phi} at which truncated density will match
         * target density
         *
         * @return the particle mass
         */
        double getParticleMass(const std::vector<std::array<double, 2>> positions,
                               const std::function<double(double, double)>
                               target_density, 
                               const std::array<double, 2> target_coords);

        template <typename DensityType>
        void initFromDensity(const DensityType& density);
        void initFromCoefficients(const std::vector<std::vector<std::complex<double>>>
                                  coefficients_);
        
        // For diagnostic purposes
        std::vector<std::vector<std::complex<double>>> getCoefficients() const;
        // Calculate norm from coefficients. This might give us an idea
        // of how well the potential is represented byt the BFE with this
        // number of terms. (eg. if the norm remains constant during N-body,
        // that's likely a good sign). Norm will be imaginary (for some reason)
        // so we return its magnitude
        double calculateAbsNorm() const;
};

// Same as method above, but object is discarded
double getParticleMass(const std::vector<std::array<double, 2>> positions,
                               const std::function<double(double, double)>
                               target_density, 
                               const std::array<double, 2> target_coords);

}