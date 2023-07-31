#include "gamma.hpp"
#include "gtest/gtest.h"

using namespace special_functions;

TEST(GammaFunctionTest, CheckIntegers) {
    ASSERT_EQ(Gamma(10), 362880);
}

TEST(GammaFunctionTest, CheckHalfIntegers) {
    ASSERT_LT(std::abs(GammaHalf(10)/1.13328e6 - 1), 1e-3);
}

// TEST(PochhammerTest, CheckIntegers) {
//     ASSERT_EQ(Pochhammer(3, 4), 360);
//     ASSERT_EQ(Pochhammer(-5, 3), -60);
// }

// TEST(PochhamerTest, CheckHalfIntegers) {
//     ASSERT_EQ(Pochhammer(3.5, 4), 563.0625);
//     ASSERT_EQ(Pochhammer(-5.5, 3), -86.625);
// }