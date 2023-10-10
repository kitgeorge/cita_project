#pragma once
#include <memory>
#include "bfe.hpp"
#include "add_functions.hpp"
#include "comp_funcs.hpp"
#include "flatten.hpp"
#include "execute_in_parallel.hpp"

namespace basis_functions {





/**
 * Calculates truncated density/potential/force functions from density
 *
 * This class uses a BFE object to find the BFE coefficients of 
 * a density object (either a function or a vector of particle positions/masses).
 * It can also accept the BFE coefficients directly. These coefficients
 * describe terms up to maximum n and |l| indices in the BFE. 
 * Density, potential and force functions are expanded in basis functions
 * and the expansion is truncated at nl_max, giving truncated density,
 * potential and force functions.
 *
 */
class PotentialFromDensity {
    const std::shared_ptr<const basis_functions::BFE> expansion;
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
        /**
         * Truncated density function of (R, phi)
         *
         * Density in (R, phi), expanded in basis functions with only
         * terms of indices n <= n_max, -l_max <= l <= l_max included.
         */
        std::function<double(double, double)> trunc_density;
        /// Truncated potential function of (R, phi)
        std::function<double(double, double)> trunc_potential;
        /// Truncated force function of (R, phi)
        std::function<std::array<double, 2>(double, double)> trunc_force;


        /**
         * Constructor
         *
         * Sets the maximum indices n and l of basis functions to be used,
         * and constructs a BFE object, a shared pointer to which is kept.
         *
         * @param n_max maximum value of n index in BFE
         * @param l_max maximum(/negative of minimum) value for l index in BFE
         *
         * @note n_max and l_max for the BFE object are set as constants in 
         * its BFETables subclass, so the values here are not connected to the 
         * BFE object. Just make sure they are not greater than the BFETables
         * values (the default values here at time of writing). This is poor
         * design but not a priority now.
         */
        PotentialFromDensity(int n_max=64, int l_max=32);
        /// Constructor as above, but stores pointer to existing BFE object
        PotentialFromDensity(const std::shared_ptr<const basis_functions::BFE> 
                             expansion_,
                             int n_max=64, int l_max=32);
        /// Copy constructor
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

        /**
         * Finds BFE coefficients -> truncated functions for a given density
         *
         * This function takes a density, either as a function or a list of particle
         * positions with masses, and uses the BFE object to find the corresponding
         * BFE coefficients (up to n_max, l_max). The truncated density, potential
         * and force functions (truncated to only the terms up to n_max, l_max) are
         * then set from these coefficients
         *
         * @param density This template is explictly instantiated for density as an
         * std::function<double(double, double)> (function of R and phi) or an
         * std::vector<std::array<double, 3>> (a vector of particles' 
         * {R, phi, m} values).
         */
        template <typename DensityType>
        void initFromDensity(const DensityType& density);
        /// Directly sets BFE coefficients -> truncated functions as above
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