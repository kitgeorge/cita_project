#include "gtest/gtest.h"
#include "theta_r_integrand.hpp"
#include "mestel.hpp"
#include "units.hpp"
#include "vector_io.hpp"

using namespace actions;

TEST(ThetaRSTQIntegrandTest, PlottingCheck) {
    double v_c = 220*Units::kms;
    double R_0 = 8*Units::kpc;
    potential::AxsymFuncs pot = potential::getMestel(v_c, R_0);
    double E = 1.05*pow(v_c, 2)/2;
    double L = R_0*v_c;
    std::function<double(double)> 
    integrand = getThetaRSTQIntegrand(pot, E, L);
    double u_max = 5;
    int N_intervals = 2e4;
    std::vector<double> integrand_values(N_intervals);
    for(int i = 0; i < N_intervals; ++i) {
        double u = -u_max + (double)i/N_intervals*2*u_max;
        integrand_values[i] = integrand(u);
    }
    utility::writeCsv("../test_data/actions/theta_r_s_t_q_integrand.csv",
                      integrand_values);
}