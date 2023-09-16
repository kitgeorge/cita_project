#include "gtest/gtest.h"
#include "sample.hpp"
#include "dehnen_df.hpp"
#include "flatten.hpp"
#include "vector_io.hpp"
#include "units.hpp"
#include "mestel.hpp"
#include "tapered_df.hpp"

using namespace df;

TEST(SampleTest, getDFSampleELPlotting) {
    double v_c = 220*Units::kms;
    double R_0 = 8*Units::kpc;
    potential::AxsymFuncs pot = potential::getMestel(v_c, R_0);
    std::function<double(double)>
    surface_density = [](double R) {
        double S_0 = 49*Units::Msun/Units::pc2;
        return S_0*exp(-(R - 8)/8);
    };
    std::function<double(double)>
    sigma_R = [](double R) {
        double sigma_R_0 = 30*Units::kms;
        return sigma_R_0*exp(-(R - 8)/8);
    };
    double E_max = pow(v_c, 2)/2 + pot.potential_R(16);
    double L_max = 16*v_c;
    double E_min = -1*E_max;
    double L_min = 0.1*L_max;
    double N_samples = 1e3;
    std::array<std::array<double, 2>, 2>
    bounds = {{ {{E_min, E_max}}, {{L_min, L_max}} }};
    std::function<double(double, double)>
    df = getNormDehnenDF(pot, surface_density, sigma_R, bounds[0]);
    std::vector<std::array<double, 2>> samples(N_samples);
    for(int i = 0; i < N_samples; ++i) {
        samples[i] = getDFSampleEL(df, bounds);
    }
    utility::writeCsv("../test_data/df/get_df_sample_e_l.csv",
                      utility::flatten(samples));

}

TEST(SampleTest, DISABLED_getDFSampleViaELDensityPlotting) {
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
    std::array<std::array<double, 2>, 2>
    E_L_bounds = {tap.E_bounds, tap.L_bounds};
    
    double u_max = 5;
    int N_u_intervals = 1000;
    int N_u_iterate = 5;

    int N_samples = 1e4;
    std::vector<std::array<double, 2>> sample_positions(N_samples);

    for(int i = 0; i < N_samples; ++i) {
        std::array<std::array<double, 2>, 2>
        sample = getDFSampleViaEL(tap.getTaperedDF(), E_L_bounds,
                                  potential::getMestel(v_c, R_0),
                                  u_max, N_u_intervals, N_u_iterate);
        sample_positions[i] = sample[0];
    }
    utility::writeCsv("../test_data/df/tapered_df_sample.csv",
                      utility::flatten(sample_positions));
}