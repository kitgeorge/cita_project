#pragma once
#include "axsym_funcs.hpp"
#include "mapXVtoAA2D.hpp"

namespace actions {

class JresFinder {
    // Variables are named differently to in Rimpei's code (eg. phi rather than psi)
    // Variables declared in order which allows constructor to work

    const potential::AxsymFuncs pot;

    const double Omega_p;

    const int N_table_intervals;

    const int N_tau;
    const int N_bisect;

    public:
        const int N_R;
        const int N_phi; 
        const double J_s_res_0;

    private:
        const double dJs;
        const std::vector<std::vector<double>> E_J;

        double calculateJsres0() const;
        std::vector<std::vector<double>> tabulateE() const;
        double getE(double J_R, double L) const;

        double calculateOmega_s(double E, double L) const;
    
    public:
        JresFinder(const potential::AxsymFuncs& pot_, int N_R_, int N_phi_,
                   double Omega_p_, int N_table_intervals_, int N_tau_,
                   int N_bisect_);
        JresFinder(const JresFinder& old);
        double getOmega_s(double J_f, double J_s) const;
        double calculateJsres(double J_f) const;
};

std::vector<std::array<double, 2>> tabulateJsres(JresFinder calc, int N_J_f_intervals);


}