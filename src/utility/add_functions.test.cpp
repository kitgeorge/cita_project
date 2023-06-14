#include "gtest/gtest.h"
#include "add_functions.hpp"

using utility::addFunctions;

TEST(AddFunctionsTest, WorksArray) {
    std::function<std::array<double, 2>(double, double, double)>
    f_0 = [] (double a, double b, double c) {
        std::array<double, 2>
        output = {a - b + c, a + b - c};
        return output;
    };
    std::function<std::array<double, 2>(double, double, double)>
    f_1 = [] (double a, double b, double c) {
        std::array<double, 2>
        output = {a + b - c, a - b + c};
        return output;
    };
    std::array<double, 2> check = addFunctions<2>({f_0, f_1}, 3.0, 4.0, 5.0);
    ASSERT_EQ(check[0], 6);
    ASSERT_EQ(check[1], 6);
}

// TEST(AddFunctionsTest, Works2F2A) {
//     std::function<double(double, double)>
//     f_0 = [=] (double a, double b) {
//         return 2*a + b;
//     };
//     std::function<double(double, double)>
//     f_1 = [=] (double a, double b) {
//         return a + 2*b;
//     };
//     ASSERT_EQ(addFunctions(f_0, f_1)(3, 4),
//               f_0(3, 4) + f_1(3, 4));
// }

// TEST(AddFunctionsTest, Works2F3A) {
//     std::function<double(double, double, double)>
//     f_0 = [=] (double a, double b, double c) {
//         return 2*a + b + 3*c;
//     };
//     std::function<double(double, double, double)>
//     f_1 = [=] (double a, double b, double c) {
//         return a + 2*b + 3*c;
//     };
//     ASSERT_EQ(addFunctions(f_0, f_1)(3, 4, 5),
//               f_0(3, 4, 5) + f_1(3, 4, 5));
// }

// TEST(AddFunctionsTest, Iterative) {
//     std::function<double(double, double)>
//     f_0 = [=] (double a, double b) {
//         return 2*a + b;
//     };
//     std::function<double(double, double)>
//     f_1 = [=] (double a, double b) {
//         return 0;
//     };
//     std::function<double(double, double)>
//     f_2 = f_0;
//     for(int i = 0; i < 100; ++i) {
//         f_2 = addFunctions(f_2, f_1);
//     }
//     ASSERT_EQ(f_2(3, 4), f_0(3, 4));
// }