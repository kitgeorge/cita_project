#include "mestel_spiral.hpp"
#include "gtest/gtest.h"
#include <optional>
#include <iostream>

using namespace potential;

class MestelSpiralTest: public ::testing::Test {
    protected:
        double v_c;
        double R_0;
        int m;
        double pattern_speed;
        double q;
        int global;
        double alpha;
        double f;

        std::optional<MestelSpiralEnvelope> envelope;
        std::optional<MestelSpiralLibrationFuncs> lib_funcs;

        void SetUp() override {
            v_c = 220*Units::kms;
            R_0 = 8*Units::kpc;
            m = 2;
            pattern_speed = v_c*(1 + sqrt(2)/m)/R_0;
            q = 0.5;
            global = 0;
            alpha = std::numbers::pi/6;
            f = 0.01;
            
            envelope.emplace(MestelSpiralEnvelope(v_c, q, global, m,
                                                  pattern_speed, alpha, f));
            lib_funcs.emplace(MestelSpiralLibrationFuncs(envelope.value()));
        }

};

TEST_F(MestelSpiralTest, EnvelopeWorks) {

    MestelSpiralEnvelope envelope(v_c, q, global, m, pattern_speed,
                                  alpha, f);
    ASSERT_EQ(envelope.s(R_0), 1);
    ASSERT_EQ(envelope.A(R_0), 0);
    ASSERT_EQ(envelope.A(R_0/(1 + sqrt(2)/2)), 1);
    ASSERT_LT(envelope.dAdR(R_0*0.9), 0);
    ASSERT_GT(envelope.dAdR(R_0*0.9/(1 + sqrt(2)/2)), 0);
    ASSERT_EQ(envelope.dAdR(R_0/(1 + sqrt(2)/2)), 0);
    ASSERT_EQ(envelope.dAdR(R_0), 0);
    ASSERT_EQ(envelope.alpha, alpha);
}

TEST_F(MestelSpiralTest, LibrationFuncsWork) {
    double J_s = 0.9*R_0*v_c/m;
    double J_f = -0.9*J_s;
    std::complex<double> clm = lib_funcs.value().clm(J_f, J_s, 1);
    double F = lib_funcs.value().F(J_f, J_s, 1);
    std::cout << F << std::endl;
}