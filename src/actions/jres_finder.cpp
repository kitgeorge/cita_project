#include "jres_finder.hpp"
#include <iostream>

namespace actions {

// Let's extend beyond J_s_res_0 in J_s, so we can also see ILR
JresFinder::JresFinder(const potential::AxsymFuncs& pot_, int N_R_, int N_phi_,
                   double Omega_p_, int N_table_intervals_, int N_tau_,
                   int N_bisect_): pot(pot_), N_bisect(N_bisect_),
                   N_R(N_R_), N_phi(N_phi_), Omega_p(Omega_p_), 
                   N_table_intervals(N_table_intervals_), N_tau(N_tau_), 
                   J_s_res_0(calculateJsres0()), dJs(1.5*J_s_res_0/N_table_intervals),
                   E_J(tabulateE()) {}


JresFinder::JresFinder(const JresFinder& old):
    pot(old.pot), N_bisect(old.N_bisect), N_R(old.N_R), N_phi(old.N_phi),
    Omega_p(old.Omega_p), N_table_intervals(old.N_table_intervals),
    N_tau(old.N_tau), J_s_res_0(old.J_s_res_0), dJs(old.dJs), 
    E_J(old.E_J) {}

double JresFinder::calculateJsres0() const {
    double R_res = pot.resonanceRadius(N_R, N_phi, Omega_p);
    return pow(R_res, 2)*pot.Omega(R_res)/N_phi;
}  

std::vector<std::vector<double>> JresFinder::tabulateE() const {
    int prog = 1;
    std::cout << N_tau << std::endl;
    std::vector<std::vector<double>> output;
    RC::gridEOverJ(&pot, N_tau, N_table_intervals + 1, N_table_intervals + 1, 
                   dJs, N_phi*dJs, output, prog);
    return output;
}

double JresFinder::getE(double J_R, double L) const {
    return RC::XGivenJ(J_R, L, dJs, N_phi*dJs, E_J);
}

double JresFinder::getOmega_s(double J_f, double J_s) const {
    double J_R = J_f + N_R*J_s;
    double L = N_phi*J_s;
    return calculateOmega_s(getE(J_R, L), L);
}

double JresFinder::calculateJsres(double J_f) const {
    std::array<double, 2> 
    bounds = {-J_f/N_R, J_s_res_0};
    assert(getOmega_s(J_f, bounds[0]) > 0);
    // assert(calculate_Omega_s(bounds[1]) < 0);
    for(int i = 0; i < N_bisect; ++i) {
        double middle = (bounds[0] + bounds[1])/2;
        if(getOmega_s(J_f, middle) > 0) {
            bounds[0] = middle;
        }
        else {
            bounds[1] = middle;
        }
    }
    return (bounds[0] + bounds[1])/2;
}

double JresFinder::calculateOmega_s(double E, double L) const {
    RC::MapXVtoAA2D calc(&pot);
    calc.mapLEtoJ(L, E, N_tau);
    return N_R*calc.wr + N_phi*(calc.wpsi - Omega_p);
}

std::vector<std::array<double, 2>> tabulateJsres(JresFinder calc, int N_J_f_intervals) {
    std::vector<std::array<double, 2>> output(N_J_f_intervals + 1);
    for(int i = 0; i < N_J_f_intervals + 1; ++i) {
        std::cout << "tabulating " << i << std::endl;
        double J_f = ((double)i/N_J_f_intervals - 1)*calc.J_s_res_0;
        output[i] = {J_f, calc.calculateJsres(J_f)};
    }
    return output;
}

}