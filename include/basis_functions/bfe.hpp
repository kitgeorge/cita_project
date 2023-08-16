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

// Basis function coefficients tabulated
// I regret tabulating this all globally. If the basis functions
// are already cached then it's not a problem, but if you have to 
// calculate them then the alpha_Ka values etc likely stay in memory,
// meaning you probably have to abort and start again (which is
// untidy even if not terribly inconvenient)

// All formulae involved (from eg Fouvry 2015) are in an anonymous namespace
// in bfe.cpp


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
        BFETables();
        BFETables(const BFETables& old);

        double getU(int n, int l, double R, double R_Ka) const;
        double getUPrime(int n, int l, double R, double R_Ka) const;
        double getD(int n, int l, double R, double R_Ka) const;
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
        BFE(double R_Ka_, int N_R_, int N_phi_);
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
        // Takes discrete density distribution (in arrays {m, R, phi})
        std::complex<double>
        getCoefficient(int n, int l,
                        const std::vector<std::array<double, 3>>& density) const;
    

};

}
