#include "gtest/gtest.h"
#include "coords.hpp"

using namespace vectors;

TEST(Coords2dTest, Consistent) {
    std::array<std::array<double, 2>, 2> cart = {{ {{1, 2}}, {{3, 4}} }};
    Coords2d point(cart, 0);
    Coords2d point2(point.polar, 1);
    for(int i = 0; i < 2; ++i) {
        for(int j = 0; j < 2; ++j) {
            ASSERT_EQ(pow(point2.cartesian[i][j] - cart[i][j], 2) < 1e-6, 1);
        }
    }
}