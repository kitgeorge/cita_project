#include "gtest/gtest.h"
#include "vector_io.hpp"

using utility::writeCsv;
using utility::readCsv;

TEST(VectorIoTest, Both1dWork) {
    std::vector<double> x = {2.3, 5.6, 2.7};
    std::string path = "../test_data/utility/both_1d_work.csv";
    writeCsv(path, x);
    std::vector<double> y = readCsv(path);
    ASSERT_EQ(x.size(), y.size());
    for(int i = 0; i < x.size(); ++i) {
        ASSERT_EQ(x[i], y[i]);
    }
}