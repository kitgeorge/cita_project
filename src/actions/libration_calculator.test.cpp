#include "gtest/gtest.h"
#include "libration_calculator.hpp"
#include "mestel_spiral.hpp"
#include "mestel.hpp"
#include "units.hpp"
#include "vector_io.hpp"
#include "flatten.hpp"
#include <csignal>
#include <memory>
#include <optional>
#include <cmath>

using namespace actions;

TEST(JacobiSpecialFunctionsTest, cnInvConsistent) {
    int N_intervals = 2000;
    JacobiSpecialFunctions jfuncs(N_intervals);
    double u = 1;
    double m = 0.5;
    ASSERT_EQ(jfuncs.cnInv(jfuncs.cn(u, m), m) - 1 < 1e-2, 1);
}

class LibrationCalculatorTest : public ::testing::Test {
    protected:
        double v_c;
        double R_0;
        double pattern_speed;
        int m;
        int l;
        double q;
        int global;
        double alpha;
        double f;

        std::optional<potential::MestelSpiralEnvelope> envelope;
        std::optional<potential::MestelSpiralLibrationFuncs> libfuncs;
        std::function<double(double, double)> F;
        std::function<double(double, double)> g;

        std::optional<JresFinder>finder;

        void SetUp() override {
            v_c = 220*Units::kms;
            R_0 = 8*Units::kpc;
            m = 2;
            pattern_speed = v_c*(1 + sqrt(2)/m)/R_0;
            q = 0.5;
            global = 0;
            alpha = std::numbers::pi/6;
            f = 0.01;
            l = 1;

            envelope.emplace(v_c, q, global, m, pattern_speed, alpha, f);
            libfuncs.emplace(envelope.value());
            // F = [this] (double J_f, double J_s_res) {
            //     return libfuncs->F(J_f, J_s_res, l);
            // };
            // g = [this] (double J_f, double J_s_res) {
            //     return libfuncs->g(J_f, J_s_res, l);
            // };
            
            potential::AxsymFuncs mestel = potential::getMestel(v_c, R_0);
            int N_table_intervals = 100;
            int N_tau = 1000;
            int N_bisect = 10;
            finder.emplace(mestel, l, m, pattern_speed, N_table_intervals, N_tau, N_bisect);
        }
};

// class LibrationCalculatorTest : public ::testing::Test {
//     protected:
//         double v_c;
//         double R_0;
//         double pattern_speed;
//         int m;
//         int l;

//         potential::MestelSpiral spiral;
//         std::function<double(double, double)> F;
//         std::function<double(double, double)> g;

//         std::unique_ptr<JresFinder> finder;

//         void SetUp() override {
//             v_c = 220*Units::kms;
//             R_0 = 8*Units::kpc;
//             pattern_speed = v_c*(1 + sqrt(2))/R_0;
//             double q = 0.5;
//             int global = 0;
//             m = 2;
//             double alpha = std::numbers::pi/6;
//             double f = 0.01;
//             l = 1;

//             spiral = potential::MestelSpiral(v_c, q, global, m,
//                                              pattern_speed, alpha, f);
//             F = [this] (double J_f, double J_s_res) {
//                 return spiral.F(J_f, J_s_res, l);
//             };
//             g = [this] (double J_f, double J_s_res) {
//                 return spiral.g(J_f, J_s_res, l);
//             };

//             potential::AxsymFuncs mestel = potential::getMestel(v_c, R_0);
//             int N_table_intervals = 100;
//             int N_tau = 1000;
//             int N_bisect = 10;
//             finder = std::make_unique<JresFinder>(mestel, l, m, pattern_speed, N_table_intervals, N_tau, N_bisect);



//         }

// };

TEST_F(LibrationCalculatorTest, DISABLED_PlottingCheck) {
    
    double J_f = -0.4*R_0*v_c;
    double delta_J_s = 0.001*R_0*v_c;
    int N_jacobi_intervals = 1000;
    double J_s_res = finder->calculateJsres(J_f);
    std::cout << J_f << ", " << J_s_res << ", " << finder->J_s_res_0 << std::endl;
    ASSERT_EQ(std::isfinite(libfuncs->F(J_f, J_s_res, l)), 1);
    ASSERT_EQ(std::isfinite(libfuncs->g(J_f, J_s_res, l)), 1);
    LibrationCalculator calc(J_f, finder.value(), libfuncs.value(), 
                             delta_J_s, N_jacobi_intervals);
    double theta_s_0 = -calc.g - std::numbers::pi/3;
    double J_s_0 = J_s_res - (J_s_res - calc.separatrices(theta_s_0)[0])/2;
    
    
    calc.setICs(J_s_0, theta_s_0);
    int N_timesteps = 1000;
    double t_max = 10*2*std::numbers::pi/calc.omega_0();
    int N_separatrix_points = 100;
    std::vector<std::array<double, 3>>
    separatrices(N_separatrix_points);
    for(int i = 0; i < N_separatrix_points; ++i) {
        double angle = ((double)i/N_separatrix_points - 0.5)*2*std::numbers::pi
                       - calc.g;
        separatrices[i][0] = angle;
        separatrices[i][1] = calc.separatrices(angle)[0];
        separatrices[i][2] = calc.separatrices(angle)[1];
        // Extremely weird bug: if we use the following for loop instead, then
        // for some reason the code continues to believe i < N_separatrix_points
        // even when i = 1200 (N_separatrix_points = 100). I don't know why 
        // getting rid of the j=2 for loop fixes it, and I find it disturbing.
        // Google test bug?
        // for(int j = 0; j < 2; ++j) {
        //     std::cout << j << std::endl;
        //     separatrices[i][j + 1] = calc.separatrices(angle)[i];
        //     std::cout << j << "a" << std::endl;
        // }
    }
    std::vector<std::array<double, 3>>
    trajectory(N_timesteps);
    for(int i = 0; i < N_timesteps; ++i) {
        double t = (double)i/N_timesteps*t_max;
        std::array<double, 2> coords = calc.trajectory(t);
        trajectory[i][0] = t;
        for(int j = 0; j < 2; ++j) {
            trajectory[i][j + 1] = coords[j];
        }
    }

    utility::writeCsv("../test_data/actions/separatrices.csv",
                      utility::flatten(separatrices));
    utility::writeCsv("../test_data/actions/libration_trajectory.csv",
                      utility::flatten(trajectory));

}