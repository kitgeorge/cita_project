
#include "gtest/gtest.h"
#include "dehnen_df.hpp"
#include "axsym_funcs.hpp"
#include "mestel.hpp"
#include "units.hpp"
#include "flatten.hpp"
#include "vector_io.hpp"

using namespace df;

TEST(DehnenDFTest, PlottingDF) {
    double v_c = 220*consts::km;
    potential::AxsymFuncs pot = potential::getMestel();
    std::function<double(double)>
    surface_density = [](double R) {
        double S_0 = 49*consts::M_solar/pow(consts::pc, 2);
        return S_0*exp(-(R - 8*consts::kpc)/(8*consts::kpc));
    };
    std::function<double(double)>
    sigma_R = [](double R) {
        double sigma_R_0 = 30*consts::km;
        return sigma_R_0*exp(-(R - 8*consts::kpc)/(8*consts::kpc));
    };
    double E_max = pow(v_c, 2)/2 + pot.potential_R(16*consts::kpc);
    double L_max = 16*consts::kpc*v_c;
    double E_min = 0.1*E_max;
    double L_min = 0.1*L_max;
    double N_E_values = 1000;
    double N_L_values = 1000;
    std::function<double(double, double)>
    df = getDehnenDF(pot, surface_density, sigma_R);
    std::vector<std::vector<double>>
    df_values(N_E_values, std::vector<double>(N_L_values));
    for(int i = 0; i < N_E_values; ++i) {
        for(int j = 0; j < N_L_values; ++j) {
            df_values[i][j] = df(E_min + (double)i/N_E_values*(E_max - E_min),
                                 L_min + (double)j/N_L_values*(L_max - L_min));
        }
    }
    utility::writeCsv("../test_data/df/dehnen_df_test.csv", 
                      utility::flatten(df_values));
}