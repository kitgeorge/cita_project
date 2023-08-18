#include "tapered_df.hpp"
#include "vector_io.hpp"
#include "flatten.hpp"
#include "gtest/gtest.h"

using namespace df;

TEST(TaperedDFTest, PlottingCheck) {
    int N_E = 1000;
    int N_L = 1000;
    std::vector<std::vector<double>>
    df_values(N_E, std::vector<double>(N_L));

    double v_c = 220*Units::kms;
    double R_0 = 1*Units::kpc;
    double active_fraction = 0.5;
    std::array<double, 2> taper_radii = {1, 10};
    std::array<double, 2> taper_indices = {4, 5};
    std::array<double, 2> cutoff_radii = {0.1, 20};
    std::function<double(double)>
    target_Q = [] (double R) {return 2;};
    TaperedDF tap(v_c, R_0, active_fraction,
                  taper_radii, taper_indices, cutoff_radii,
                  target_Q);

    double delta_E = (tap.E_bounds[1] - tap.E_bounds[0])/N_E;
    double delta_L = (tap.L_bounds[1] - tap.L_bounds[0])/N_L;
    for(int i = 0; i < N_E; ++i) {
        double E = tap.E_bounds[0] + i*delta_E;
        for(int j = 0; j < N_L; ++j) {
            double L = tap.L_bounds[0] + j*delta_L;
            df_values[i][j] = tap.getTaperedDF()(E, L);
        }
    }
    utility::writeCsv("../test_data/df/tapered_df.csv",
                      utility::flatten(df_values));
}