#include "gtest/gtest.h"
#include "bfe_pochhammer.hpp"
#include <functional>
#include <cmath>

using namespace special_functions;

TEST(PochhammerTest, Consistent) {
    ASSERT_EQ(Pochhammer(3, 4), getPochhammerInt(3, 4));
    ASSERT_EQ(Pochhammer(-3, 4), getPochhammerInt(-3, 4));
    ASSERT_EQ(Pochhammer(3.5, 4), getPochhammerHalfInt(3, 4));
}