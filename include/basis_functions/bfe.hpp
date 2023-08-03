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

namespace {
static constexpr int k_max = 10;
static constexpr int l_max = 50;
static constexpr int n_max = 50;
static constexpr int i_max = 10;
static constexpr int j_max = 50;

static constexpr int N_R_tabulated = 1e4;
}

utility::vector5d<LooongDouble> getAlphaKaValues();
utility::vector4d<LooongDouble> getBetaKaValues();
utility::vector3d<double> getPValues();
utility::vector3d<double> getSValues();

utility::vector5d<LooongDouble>& 
alpha_Ka_values();
utility::vector4d<LooongDouble>& 
beta_Ka_values();
utility::vector3d<double>& 
P_values();
utility::vector3d<double>& 
S_values();

LooongDouble getAlphaKa(int k, int l, int n, int i, int j);
LooongDouble getBetaKa(int k, int l, int n, int j);
double getP(int k, int l, int n);
double getS(int k, int l, int n);


utility::vector4d<double> getUValues();
utility::vector4d<double> getUPrimeValues();
utility::vector4d<double> getDValues();

utility::vector4d<double> calculateUValues();
utility::vector4d<double> calculateUPrimeValues();
utility::vector4d<double> calculateDValues();

utility::vector4d<double>& U_values();
utility::vector4d<double>& UPrime_values();
utility::vector4d<double>& D_values();

double getU(int k, int n, int l, double R, double R_Ka);
double getUPrime(int k, int n, int l, double R, double R_Ka);
double getD(int k, int n, int l, double R, double R_Ka);

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
    const int k_Ka;
    const double R_Ka;
    const int N_R;
    const int N_phi;

    // double U(int n, int l, double R) const;
    // double UPrime(int n, int l, double R) const;
    // double D(int n, int l, double R) const;

    public:
        BFE(int k_Ka_, double R_Ka_, int N_R_, int N_phi_);
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
