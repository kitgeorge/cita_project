#include "gtest/gtest.h"
#include "jres_finder.hpp"
#include "mestel.hpp"
#include "vector_io.hpp"
#include "flatten.hpp"

using namespace actions;

TEST(JresFinderTest, DISABLED_PlottingTabulated) {
    double v_c = 220*Units::kms;
    double R_0 = 8*Units::kpc;
    potential::AxsymFuncs pot = potential::getMestel(v_c, R_0);
    int N_R = 1;
    int N_phi = 2;
    double Omega_p = pot.Omega(R_0) + pot.kappa(R_0)/N_phi;
    int N_table_intervals = 100;
    int N_tau = 1000;
    int N_bisect = 10;
    int N_J_f_intervals = 100;
    JresFinder calc(pot, N_R, N_phi, Omega_p, N_table_intervals,
                    N_tau, N_bisect);
    std::vector<std::array<double, 2>>
    Jres_values = tabulateJsres(calc, N_J_f_intervals);
    std::vector<std::array<double, 2>>
    L_J_R_res(N_J_f_intervals + 1);
    for(int i = 0; i < N_J_f_intervals + 1; ++i) {
        L_J_R_res[i][0] = N_phi*Jres_values[i][1];
        L_J_R_res[i][1] = N_R*Jres_values[i][1] + Jres_values[i][0];
    }
    utility::writeCsv("../test_data/actions/jres_finder.csv",
                      utility::flatten(L_J_R_res));
}