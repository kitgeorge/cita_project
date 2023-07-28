#pragma once
#include "gamma.hpp"
#include "units.hpp"
#include "comp_funcs.hpp"
#include <complex>
#include <cassert>


namespace basis_functions {

std::complex<double>
scalarProduct(const std::function<std::complex<double>(double, double)>& pot_conj,
               const std::function<std::complex<double>(double, double)>& density,
               double R_max, int N_R, int N_phi);

std::complex<double>
scalarProduct(const std::function<std::complex<double>(double, double)>& pot_conj,
               const std::vector<std::array<double, 3>>& density, double R_max);

class BFE {
    // Following notation in Fouvry et al (2015)
    const int k_Ka;
    const double R_Ka;
    const int N_R;
    const int N_phi;

    double alpha_Ka(int k, int l, int n, int i, int j) const;
    double beta_Ka(int k, int l, int n, int j) const;
    double P(int k, int l, int n) const;
    double S(int k, int l, int n) const;
    double U(int n, int l, double R) const;
    double UPrime(int n, int l, double R) const;
    double D(int n, int l, double R) const;

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
