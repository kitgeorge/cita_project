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

TEST(NDVectorsTest, ReshapeWorks4d) {
    std::array<int, 4> shape = {2, 3, 4, 5};
    int N = shape[0]*shape[1]*shape[2]*shape[3];
    std::vector<int> flat(N);
    for(int i = 0; i < N; ++i) {
        flat[i] = 2*i;
    }
    vector4d<int> result = reshape(flat, shape);
    for(int i = 0; i < shape[0]; ++i) {
        for(int j = 0; j < shape[1]; ++j) {
            for(int k = 0; k < shape[2]; ++k) {
                for(int l = 0; l < shape[3]; ++l) {
                    ASSERT_EQ(result[i][j][k][l], 
                              flat[i*shape[1]*shape[2]*shape[3]
                                   + j*shape[2]*shape[3] 
                                   + k*shape[3] + l]);
                }
            }
        }
    }
}