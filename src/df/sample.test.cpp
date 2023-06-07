#include "gtest/gtest.h"
#include "sample.hpp"
#include "dehnen_df.hpp"
#include "flatten.hpp"
#include "vector_io.hpp"
#include "units.hpp"
#include "mestel.hpp"

using namespace df;

TEST(SampleTest, getDFSampleELPlotting) {
    double v_c = 220*Units::kms;
    potential::AxsymFuncs pot = potential::getMestel(v_c);
    std::function<double(double)>
    surface_density = [](double R) {
        double S_0 = 49*Units::M_sun/Units::pc2;
        return S_0*exp(-(R - 8)/8);
    };
    std::function<double(double)>
    sigma_R = [](double R) {
        double sigma_R_0 = 30*Units::kms;
        return sigma_R_0*exp(-(R - 8)/8);
    };
    double E_max = pow(v_c, 2)/2 + pot.potential_R(16);
    double L_max = 16v_c;
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