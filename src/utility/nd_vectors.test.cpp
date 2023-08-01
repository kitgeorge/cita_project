#include "nd_vectors.hpp"
#include "gtest/gtest.h"

using namespace utility;

TEST(NDVectorsTest, Works5d) {
    std::array<int, 5> shape = {2, 3, 4, 5, 6};
    vector5d<double> vec = makeShape<double>(shape);
    ASSERT_EQ(vec.size(), shape[0]);
    ASSERT_EQ(vec[0].size(), shape[1]);
    ASSERT_EQ(vec[0][0].size(), shape[2]);
    ASSERT_EQ(vec[0][0][0].size(), shape[3]);
    ASSERT_EQ(vec[0][0][0][0].size(), shape[4]);
}