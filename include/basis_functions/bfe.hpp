#pragma once
#include "gamma.hpp"
#include "bfe_pochhammer.hpp"
#include "units.hpp"
#include "comp_funcs.hpp"
#include "nd_vectors.hpp"
#include "vector_io.hpp"
#include "flatten.hpp"
#include "execute_in_parallel.hpp"
#include <complex>
#include <cassert>
#include <boost/multiprecision/gmp.hpp>
// With such high precision, I'll probably have to run on server for RAM

using LooongDouble = boost::multiprecision::mpf_float_1000;

namespace basis_functions {

// All formulae involved (from eg Fouvry 2015) are in an anonymous namespace
// in bfe.cpp

/**
 * Tabulates basis functions in R
 *
 * R-dependent parts of the Kalnajs basis functions for a 2D disc 
 * (notation here follows Fouvry 2015) are tabulated. They are either
 * read from a text file, or calculated (which involves temporary 
 * tabulation of other functions) and cached in said text file for 
 * future use.
 *
 * @note maximum values of indices n, l and sum indices i, j, as
 * well as the value of basis label k_Ka and the number of R
 * values to tabulate over, are set in the class declaration as
 * constants. eg. 0 <= n <= n_max. l is non-negative here, but
 * negative-l coefficients are simply the complex conjugates of
 * their positive-l counterparts.
 *
 * @note values cached in text files in "../cache/basis_functions/"
 *
 * @note The alpha_Ka values and beta_Ka values are used as
 * coefficients of very high-order polynomials. I've so far taken
 * a brute-force approach, calculating these values to high precision
 * as 1000-digit LooongDoubles. However, if the instability to 
 * precision here arises only from the high-order polynomials,
 * we could perhaps get around this by putting the polynomials in 
 * Lagrange form. However, the alpha_Ka and beta_Ka values themselves
 * may be dramatically varying in magnitude and so require high precision
 * to cancel correctly (as they are functions of Pochhammer symbols).
 * I should probably investigate at some point to see if we can get
 * rid of the brute-force approach.
 */
class BFETables {
    static constexpr int k_Ka = 10;
    static constexpr int l_max = 32;
    static constexpr int n_max = 64;
    static constexpr int i_max = 10;
    static constexpr int j_max = 64;

    static constexpr int N_R_tabulated = 1e4;

    std::optional<special_functions::GammaTables> gtables;
    std::optional<special_functions::PochhammerTables> ptables;

    std::optional<utility::vector4d<LooongDouble>> alpha_Ka_values;
    std::optional<utility::vector3d<LooongDouble>> beta_Ka_values;
    std::optional<utility::vector2d<double>> P_values;
    std::optional<utility::vector2d<double>> S_values;

    const utility::vector3d<double> U_values;
    const utility::vector3d<double> UPrime_values;
    const utility::vector3d<double> D_values;

    void ensureSubTablesExist();
    utility::vector4d<LooongDouble> getAlphaKaValues() const;
    utility::vector3d<LooongDouble> getBetaKaValues() const;
    utility::vector2d<double> 
    getPSValues(std::function<double(int, int, int)> which) const;
    utility::vector2d<double> getPValues() const;
    utility::vector2d<double> getSValues() const;

    double P(int k, int l, int n) const;
    double S(int k, int l, int n) const;
    LooongDouble alpha_Ka(int k, int l, int n, int i, int j) const;
    LooongDouble beta_Ka(int k, int l, int n, int j) const;

    LooongDouble getAlphaKa(int l, int n, int i, int j) const;
    LooongDouble getBetaKa(int l, int n, int j) const;
    double getP(int l, int n) const;
    double getS(int l, int n) const;

    std::optional<utility::vector3d<double>> 
    readUUprimeDValues(std::string path);
    utility::vector3d<double> getUValues();
    utility::vector3d<double> getUPrimeValues();
    utility::vector3d<double> getDValues();

    utility::vector3d<double> 
    calculateUUpDValues( 
            std::function<double(int, int, int, double)> which);
    utility::vector3d<double> calculateUValues();
    utility::vector3d<double> calculateUPrimeValues();
    utility::vector3d<double> calculateDValues();

    double U(int k, int n, int l, double R_norm) const;
    double UPrime(int k, int n, int l, double R_norm) const;
    double D(int k, int n, int l, double R_norm) const;


    public:
        /**
         * Tabulates functions
         *
         * Checks "../cache/basis_functions/" for cached basis function
         * files (and checks that each is of the right length for our
         * chosen parameters). If so, it reads them and tabulates. 
         * Otherwise, it calculates and tabulates, and writes to these
         * files (overwriting if applicable).
         *
         * @note intermediate functions are tabulated in std::optional
         * member variables, which are afterwards cleared to free up
         * memory (in case usage was excessive).
         */
        BFETables();
        /**
         * Copy constructor
         *
         * Copies tabulated functions (but leaves alone std::optional
         * tables). 
         *
         * @param old copied object
         */
        BFETables(const BFETables& old);

        /**
         * Retrieves a U function value from table for a given R
         *
         * @param n n-index of U basis function
         * @param l l-index of U basis function
         * @param R value of R to retrieve value for
         * @param R_Ka R_Ka value chosen for basis function set
         * (perhaps this last one should be refactored as a member variable
         * as it helps define all the basis functions)
         *
         * @note U values are tabulated over N_R_tabulated bins in R,
         * between 0 and R_Ka. The central U value of each bin is
         * tabulated. This function retrieves that central U value
         * for the bin which contains R (there is no interpolation).
         */
        double getU(int n, int l, double R, double R_Ka) const;
        /// Retrieves a U_prime value, as getU does U
        double getUPrime(int n, int l, double R, double R_Ka) const;
        /// Retrieves a D value, as getU does U
        double getD(int n, int l, double R, double R_Ka) const;

        /**
         * Gives U_{n, l} function for given n, l
         *
         * @param n n-index of U basis function
         * @param l l-index of U basis Function
         * @return an std::function which contains tabuled
         * values of U_{n, l}, only for the speicified n and l,
         * takes as inputs R and scale radius R_Ka of the basis
         * functions, and returns a value for U_{n, l}.
         *
         * @note Doing it this way allows us not to copy the entire
         * BFETables object when packaging the full basis function
         * into an std::function.
         *
         * @note As above, the U_{n, l} values are tabulated over N_R_tabulated
         * bins in R, between 0 and R_Ka. The central U_{n, l} value of each
         * bin is tabulated, and the value corresponding to the bin which
         * contains R is what is returned by the std::function.
         */
        std::function<double(double, double)> 
        getUFunction(int n, int l) const;
        // As above but for UPrime
        std::function<double(double, double)> 
        getUPrimeFunction(int n, int l) const;
        // As above but for D
        std::function<double(double, double)> 
        getDFunction(int n, int l) const;
};

////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////

std::complex<double>
scalarProduct(const std::function<std::complex<double>(double, double)>& pot_conj,
               const std::function<std::complex<double>(double, double)>& density,
               double R_max, int N_R, int N_phi);

std::complex<double>
scalarProduct(const std::function<std::complex<double>(double, double)>& pot_conj,
               const std::vector<std::array<double, 3>>& density, double R_max);

////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////


class BFE {
    // Following notation in Fouvry et al (2015)
    const double R_Ka;
    // Integration numbers for scalar product (continuous density)
    const int N_R;
    const int N_phi;

    const BFETables tables;

    // double U(int n, int l, double R) const;
    // double UPrime(int n, int l, double R) const;
    // double D(int n, int l, double R) const;

    public:
        BFE(double R_Ka_=20*Units::kpc, int N_R_=5000, int N_phi_=5000);
        // BFE(double R_Ka_=20*Units::kpc, int N_R_=1000, int N_phi_=1000);
        BFE(const BFE& old);

        // Basis functions
        // I'll try these with a lambda capture of *this, which is a bit
        // weird to think about but hopefully should fully wrap up the
        // basis functions in these std::function containers
        std::function<std::complex<double>(double, double)> psi(int n, int l) const;
        std::function<std::complex<double>(double, double)> rho(int n, int l) const;
        std::function<std::array<std::complex<double>, 2>(double, double)> 
        psi_f(int n, int l) const;

        std::complex<double> 
        getCoefficient(int n, int l, 
                        const std::function<double(double, double)>& 
                        density) const;
        // Takes discrete density distribution (in arrays {R, phi, m})
        std::complex<double>
        getCoefficient(int n, int l,
                        const std::vector<std::array<double, 3>>& density) const;
    

};

}
