#include "gtest/gtest.h"
#include "execute_in_parallel.hpp"

TEST(ExecuteInParallelTest, IntWorks) {
    int N_functions = 1000;
    std::vector<std::function<int()>> functions(N_functions);
    for(int i = 0; i < N_functions; ++i) {
        functions[i] = [i]() {
            return i*i;
        };
    }
    std::vector<int> results = multithreading::executeInParallel(functions);
    ASSERT_EQ(results.size(), N_functions);
    for(int i = 0; i < N_functions; ++i) {
        ASSERT_EQ(results[i], i*i);
    }
}