#include "gtest/gtest.h"
#include "flatten.hpp"

using utility::flatten;

TEST(FlattenTest, DoubleWorks2d) {
    std::vector<std::vector<double>> x = {{1.3, 2.6, 3.2}, {7.6, 5.7, 3.8}};
    std::vector<double> y = flatten(x);
    ASSERT_EQ(y.size(), x.size()*x[0].size());
    for(int i = 0; i < x.size(); ++i) {
        for(int j = 0; j < x[0].size(); ++j) {
            ASSERT_EQ(x[i][j], y[i*x[0].size() + j]);
        }
    }
}

TEST(FlattenTest, IntWorks2d) {
    std::vector<std::vector<int>> x = {{2, 6, 3}, {8, 2, 6}};
    std::vector<int> y = flatten(x);
    ASSERT_EQ(y.size(), x.size()*x[0].size());
    for(int i = 0; i < x.size(); ++i) {
        for(int j = 0; j < x[0].size(); ++j) {
            ASSERT_EQ(x[i][j], y[i*x[0].size() + j]);
        }
    }
}

TEST(FlattenTest, DoubleWorks3d) {
    std::vector<std::vector<std::vector<double>>>
    x = {{{3.5, 2.7, 4.2, 4.8}, {2.7, 3.2, 3.6, 1.2}, {2.6, 7.8, 4.3, 2.6}},
         {{2.5, 7.8, 5.4, 3.2}, {2.6, 8.5, 7.3, 2.6}, {7.4, 3.6, 2.4, 6.3}}};
    std::vector<double> y = flatten(x);
    ASSERT_EQ(y.size(), x.size()*x[0].size()*x[0][0].size());
    for(int i = 0; i < x.size(); ++i) {
        for(int j = 0; j < x[0].size(); ++j) {
            for(int k = 0; k < x[0][0].size(); ++k) {
                ASSERT_EQ(x[i][j][k],
                          y[i*x[0].size()*x[0][0].size() + j*x[0][0].size() + k]);
            }
        }
    }
}

TEST(FlattenTest, IntWorks4d) {
    std::vector<std::vector<std::vector<std::vector<int>>>> 
    x = {{{{3, 5, 6, 2}, {3, 6, 2, 7}, {7, 3, 2, 8}},
          {{4, 7, 3, 8}, {3, 8, 1, 3}, {6, 8, 4, 2}}}};
    std::vector<int> y = flatten(x);
    std::array<int, 4> shape = utility::getShape(x);
    ASSERT_EQ(y.size(), shape[0]*shape[1]*shape[2]*shape[3]);
    for(int i = 0; i < shape[0]; ++i) {
        for(int j = 0; j < shape[1]; ++j) {
            for(int k = 0; k < shape[2]; ++k) {
                for(int l = 0; l < shape[3]; ++l) {
                    ASSERT_EQ(x[i][j][k][l], y[i*shape[1]*shape[2]*shape[3]
                                               + j*shape[2]*shape[3]
                                               + k*shape[3] + l]);
                }
            }
        }
    }

}

