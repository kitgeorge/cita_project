#include "gtest/gtest.h"
#include "add_functions.hpp"

using utility::addFunctions;

TEST(AddFunctionsTest, Works2F2A) {
    std::function<double(double, double)>
    f_0 = [=] (double a, double b) {
        return 2*a + b;
    };
    std::function<double(double, double)>
    f_1 = [=] (double a, double b) {
        return a + 2*b;
    };
    ASSERT_EQ(addFunctions(f_0, f_1)(3, 4),
              f_0(3, 4) + f_1(3, 4));
}

TEST(AddFunctionsTest, Works2F3A) {
    std::function<double(double, double, double)>
    f_0 = [=] (double a, double b, double c) {
        return 2*a + b + 3*c;
    };
    std::function<double(double, double, double)>
    f_1 = [=] (double a, double b, double c) {
        return a + 2*b + 3*c;
    };
    ASSERT_EQ(addFunctions(f_0, f_1)(3, 4, 5),
              f_0(3, 4, 5) + f_1(3, 4, 5));
}