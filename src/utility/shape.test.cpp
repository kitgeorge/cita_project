#include "gtest/gtest.h"
#include "shape.hpp"

TEST(MakeShapeTest, Works3d) {
    std::array<int, 3> shape = {2, 3, 4};
    std::vector<std::vector<std::vector<int>>>
    hello = utility::makeShape<int>(shape);
    ASSERT_EQ(hello.size(), 2);
    for(int i = 0; i < 2; ++i) {
        ASSERT_EQ(hello[i].size(), 3);
        for(int j = 0; j < 3; ++j) {
            ASSERT_EQ(hello[i][j].size(), 4);
        }
    }
}

TEST(MakeShapeTest, Works4d) {
    std::array<int, 4> shape = {2, 3, 4, 5};
    std::vector<std::vector<std::vector<std::vector<int>>>>
    hello = utility::makeShape<int>(shape);
    ASSERT_EQ(hello.size(), 2);
    for(int i = 0; i < 2; ++i) {
        ASSERT_EQ(hello[i].size(), 3);
        for(int j = 0; j < 3; ++j) {
            ASSERT_EQ(hello[i][j].size(), 4);
            for(int k = 0; k < 4; ++k) {
                ASSERT_EQ(hello[i][j][k].size(), 5);
            }
        }
    }
}

TEST(GetShapeTest, Consistent3d) {
    std::array<int, 3> shape = {2, 3, 4};
    std::vector<std::vector<std::vector<int>>> hello = utility::makeShape<int>(shape);
    std::array<int, 3> shape_ = utility::getShape(hello);
    for(int i = 0; i < 3; ++i) {
        ASSERT_EQ(shape[i], shape_[i]);
    }
}

TEST(GetShapeTest, Consistent4d) {
    std::array<int, 4> shape = {2, 3, 4, 5};
    std::vector<std::vector<std::vector<std::vector<int>>>> 
    hello = utility::makeShape<int>(shape);
    std::array<int, 4> shape_ = utility::getShape(hello);
    for(int i = 0; i < 4; ++i) {
        ASSERT_EQ(shape[i], shape_[i]);
    }
}
