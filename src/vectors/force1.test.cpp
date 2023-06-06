#include "gtest/gtest.h"
#include "force1.hpp"

using vectors::Force;

class ForceTest : public ::testing::Test {
    protected:
        void SetUp() override {
            comps_0 = {3.4, 3.7, 1.6};
            comps_1 = {4.6, 2.5, 3.9};
            f_0.set_components(comps_0);
            f_1.set_components(comps_1);
        }
        std::array<double, 3> comps_0;
        std::array<double, 3> comps_1;
        Force f_0;
        Force f_1;
};

TEST_F(ForceTest, SetWorks) {
    ASSERT_EQ(f_0.f_R, comps_0[0]);
    ASSERT_EQ(f_0.f_phi, comps_0[1]);
    ASSERT_EQ(f_0.f_z, comps_0[2]);
    ASSERT_EQ(f_1.f_R, comps_1[0]);
    ASSERT_EQ(f_1.f_phi, comps_1[1]);
    ASSERT_EQ(f_1.f_z, comps_1[2]);    
}

TEST_F(ForceTest, MultiplyWorks) {
    double scalar = 2.7;
    f_0.multiply(scalar);
    for(int i = 0; i < 3; ++i) {
        comps_0[i] = scalar*comps_0[i];
    }
    ASSERT_EQ(f_0.f_R, comps_0[0]);
    ASSERT_EQ(f_0.f_phi, comps_0[1]);
    ASSERT_EQ(f_0.f_z, comps_0[2]);
}

TEST_F(ForceTest, AddWorks) {
    f_0.add(f_1);
    for(int i = 0; i < 3; ++i) {
        comps_0[i] += comps_1[i];
    }
    ASSERT_EQ(f_0.f_R, comps_0[0]);
    ASSERT_EQ(f_0.f_phi, comps_0[1]);
    ASSERT_EQ(f_0.f_z, comps_0[2]);
}

TEST_F(ForceTest, OperatorsWork) {
    Force f_3 = f_1;
    f_3.add(f_0);
    Force f__3 = f_1 + f_0;
    Force f___3 = f_1;
    f___3 += f_0;
    ASSERT_EQ(f_3.f_R, f__3.f_R);
    ASSERT_EQ(f_3.f_phi, f__3.f_phi);
    ASSERT_EQ(f_3.f_z, f__3.f_z);
    ASSERT_EQ(f_3.f_R, f___3.f_R);
    ASSERT_EQ(f_3.f_phi, f___3.f_phi);
    ASSERT_EQ(f_3.f_z, f___3.f_z);

    Force f_4 = f_1;
    Force f_5 = f_0;
    f_5.multiply(-1);
    f_4.add(f_5);
    Force f__4 = f_1;
    f__4 -= f_0;
    ASSERT_EQ(f_4.f_R, f__4.f_R);
    ASSERT_EQ(f_4.f_phi, f__4.f_phi);
    ASSERT_EQ(f_4.f_z, f__4.f_z);

    double scalar = 4.7;
    Force f_6 = f_1;
    f_6.multiply(scalar);
    Force f__6 = f_1;
    f__6 *= scalar;
    Force f_7 = f_1;
    f_7.multiply(1/scalar);
    Force f__7 = f_1;
    f__7 /= scalar;
    ASSERT_EQ(f_6.f_R, f__6.f_R);
    ASSERT_EQ(f_6.f_phi, f__6.f_phi);
    ASSERT_EQ(f_6.f_z, f__6.f_z);
    ASSERT_EQ(f_7.f_R, f__7.f_R);
    ASSERT_EQ(f_7.f_phi, f__7.f_phi);
    ASSERT_EQ(f_7.f_z, f__7.f_z);
}