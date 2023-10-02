#include "gtest/gtest.h"
#include "theta_r_integrator.hpp"
#include "units.hpp"
#include "mestel.hpp"
#include "vector_io.hpp"
#include "flatten.hpp"
#include <iostream>
#include <memory>

using namespace actions;

class ThetaRIntegratorTest: public ::testing::Test {
    protected:
        double E;
        double L;
        double R_0;
        std::unique_ptr<ThetaRIntegrator> integrator;
        std::unique_ptr<potential::AxsymFuncs> pot;

        void SetUp() override {
            double v_c = 220*Units::kms;
            R_0 = 8*Units::kpc;
            pot = std::make_unique<potential::AxsymFuncs>(potential::getMestel(v_c, R_0));
            L = v_c*R_0;
            E = 1.01*pow(v_c, 2)/2;
            double u_max = 5;
            int N_intervals = 1e3;
            int N_iterate = 10;
            integrator = std::make_unique<ThetaRIntegrator>(*pot, E, L,
                                                        u_max, N_intervals,
                                                        N_iterate);
        }
};

TEST_F(ThetaRIntegratorTest, TVsKappa) {
    double kappa = pot->kappa(R_0);
    ASSERT_EQ(pow(integrator->getT()/(2*std::numbers::pi/kappa) - 1, 2) < 1e-4, 1);
}

TEST_F(ThetaRIntegratorTest, CoordsPlotting) {
    int N_angles = 100;
    std::vector<std::array<double, 2>> coords(N_angles);
    for(int i = 0; i < N_angles; ++i) {
        double angle = (double)i/N_angles*2*std::numbers::pi;
        coords[i] = integrator->getCoords(angle);
        coords[i][1] = coords[i][1]*Units::kms_i;
    }
    utility::writeCsv("../test_data/actions/theta_r_coords.csv",
                      utility::flatten(coords));
}

TEST_F(ThetaRIntegratorTest, CoordsConsistent) {
    int N_angles = 100;
    std::vector<double> angle_errors;
    for(int i = 0; i < N_angles; ++i) {
        double angle = (double)i/N_angles*2*std::numbers::pi;
        angle_errors[i] = integrator->calculateThetaR(integator->getCoords(angle)) - angle;
        EXPECT_LT(angle_errors[i], 1e-2);
    }
}