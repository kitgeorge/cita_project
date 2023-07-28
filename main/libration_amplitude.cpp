#include "mapXVtoAA2D.hpp"
#include "axsym_funcs.hpp"
#include "vector_io.hpp"
#include "mestel.hpp"
#include <numbers>
#include <iostream>


int main() {
    double v_c = 220*Units::kms;
    double R_0 = 8*Units::kpc;
    double L_0 = R_0*v_c;
    double delta_L = 1e-5*L_0;
    int N_J_R_values = 100;
    std::vector<double> J_R_values(N_J_R_values);
    for(int i = 0; i < N_J_R_values; ++i) {
        J_R_values[i] = (double)i/N_J_R_values*L_0/10;
    }
    std::vector<double> derivative_values(N_J_R_values);

    potential::AxsymFuncs pot = potential::getMestel(v_c, R_0);
    
    std::vector<std::vector<double>> E_J;
    bool prog = 1;
    int N_tau = 1e3;
    RC::gridEOverJ(&pot, N_tau, 1002, 5, L_0/1e4, delta_L,
                   E_J, prog, L_0 - 2*delta_L);
 

    RC::MapXVtoAA2D calc(&pot);

    for(int i = 0; i < N_J_R_values; ++i) {
        std::cout << i << std::endl;
        std::array<double, 2> E_values;
        E_values[0] = RC::XGivenJ(J_R_values[i] - delta_L, L_0 - delta_L,
                                  L_0/1e4, delta_L, E_J, L_0 - 2*delta_L);
        E_values[1] = RC::XGivenJ(J_R_values[i] + delta_L, L_0 + delta_L,
                                  L_0/1e4, delta_L, E_J, L_0 - 2*delta_L);
        std::cout << E_values[0] << ", " << E_values[1] << std::endl;
        std::array<double, 2> Omegas;
        calc.mapLEtoJ(L_0 - delta_L, E_values[0], N_tau);
        Omegas[0] = 2*std::numbers::pi*(1/calc.Tr + 1/calc.Tpsi);
        calc.mapLEtoJ(L_0 + delta_L, E_values[1], N_tau);
        Omegas[1] = 2*std::numbers::pi*(1/calc.Tr + 1/calc.Tpsi);
        derivative_values[i] = (Omegas[1] - Omegas[0])/(2*delta_L);
    }
    utility::writeCsv("../data/libration_amplitude.csv", derivative_values);
}