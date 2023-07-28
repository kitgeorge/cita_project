#pragma once
#include "units.hpp"
#include "mestel.hpp"
#include "jres_finder.hpp"
#include "execute_in_parallel.hpp"
#include <complex>
#include <array>
#include <numbers>
#include <optional>


// This class takes a Mestel v_c and disc mass fraction q.
// It generates a spiral potential for a given m number,
// pattern speed, pitch angle alpha and fraction of the disc's surface density f.
// A PotentialFuncs object is created for use in test particle
// simulations. F and g std::functions of (J_f, J_s_res, l) are created for libration calculations.
// All these can be either for a "global" spiral with no additional
// envelope (same surface density fraction everywhere) or
// a local spiral with an envelope cos(s) inside the Lindblad resonances and
// 0 elsewhere. It'll be interesting to see the difference in 
// heating between the two patterns

namespace potential {

class MestelSpiralEnvelope {

    public:
        const double v_c;
        const double R_0;
        const double Omega_p;
        const int m;
        const double alpha;
        const double f;
        const int global;
        const double Sigma_0;

        MestelSpiralEnvelope(double v_c_, double q, int global_,
                             int m_, double Omega_p_, double alpha_,
                             double f_);    
        MestelSpiralEnvelope(const MestelSpiralEnvelope& old);

        double s(double R) const;
        double A(double R) const;
        double dAdR(double R) const;

};

// Obtains a density to which the potential we might get from other 
// functions here does NOT correspond. They share the same A(R) envelope
// for convenience so we can obtain either a vaguely realisitic spiral
// density or potential
std::function<double(double, double, double)>
getMestelSpiralDensity(const MestelSpiralEnvelope& envelope);

class LibrationFuncs {
    public:
        virtual std::complex<double> 
        clm(double J_f, double J_s, int l) const=0;
        virtual double 
        F(double J_f, double J_s, int l) const=0;
        virtual double 
        g(double J_f, double J_s, int l) const=0;
};

class MestelSpiralLibrationFuncs : public LibrationFuncs {

    const MestelSpiralEnvelope envelope;

    double Rg(double J_phi) const;
    double Omega(double R) const;
    double kappa(double R) const;
    double e(double J_R, double J_phi) const;
    double gamma(double J_phi) const;


    public:
        MestelSpiralLibrationFuncs(const MestelSpiralEnvelope& envelope_);
        MestelSpiralLibrationFuncs(const MestelSpiralLibrationFuncs& old);


        std::complex<double> clm(double J_f, double J_s, int l) const;
        double F(double J_f, double J_s_res, int l) const;
        double g(double J_f, double J_s_res, int l) const;
      

};

class MestelSpiralPotential {
    const MestelSpiralEnvelope envelope;
    public:
        MestelSpiralPotential(const MestelSpiralEnvelope& envelope_);
        MestelSpiralPotential(const MestelSpiralPotential& old);
        double potential(double R, double phi, double t) const;
        std::array<double, 2> force(double R, double phi, double t) const;

};

struct MestelObjects {
    const int N_phi;
    const AxsymFuncs mestel;
    const MestelSpiralEnvelope envelope;
    const MestelSpiralLibrationFuncs libfuncs;
    const MestelSpiralPotential spiral;
    const std::array<double, 3> R_res;
    const std::array<std::optional<actions::JresFinder>, 3>  
    res_finders;
    MestelObjects(double v_c, double R_0, double Omega_p, double q,
                  int global, double alpha, double f, int N_phi_,
                  int N_table_intervals, int N_tau, int N_bisect);
    MestelObjects(const MestelSpiralEnvelope& envelope_, 
                  int N_table_intervals, int N_tau, int N_bisect);
    MestelObjects(const MestelObjects& old);
    private:
        std::array<std::optional<actions::JresFinder>, 3> 
        constructFinders(int N_phi, double Omega_p,
                         int N_table_intervals, int N_tau,
                         int N_bisect) const;
};

MestelSpiralEnvelope getStandardEnvelope(int global);

std::function<double(double, double, double)> getStandardSpiralDensity();
std::function<double(double, double, double)> getGlobalSpiralDensity();
MestelObjects getStandardSpiral();
MestelObjects getGlobalSpiral();

}