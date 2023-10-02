#include "simple_timer.hpp"
#include "gtest/gtest.h"
#include <cstdlib>

using namespace utility;

TEST(SimpleTimerTest, DISABLED_CheckDoubleReadTime) {
    SimpleTimer timer;
    // timer.start();
    // timer.stop();
    // std::cout << "Timer background: " << timer.getDuration_ns()
    //           << "ns" << std::endl;
    // double x = 34543.2523543e30;
    // timer.start();
    // int y = 3.532;
    // timer.stop();
    // std::cout << "Init y: " << timer.getDuration_ns() << "ns" << std::endl;
    // timer.start();
    // y = x;
    // timer.stop();
    // std::cout << "Read double: " << timer.getDuration_ns() << "ns" << std::endl;
    const int N_a = 1000000;
    std::vector<double> a(N_a);
    for(int i = 0; i < N_a; ++i) {
        a[i] = (double)std::rand()/RAND_MAX;
    }
    std::vector<double> b(N_a);
    timer.start();
    for(int i = 0; i < N_a; ++i) {
        b[i] = a[i];
    }
    timer.stop();
    for(int i = 0; i < N_a; ++i) {
        std::cout << b[i];
    }
    std::cout << std::endl << "Read double: " << timer.getDuration_ns()/N_a
              << "ns" << std::endl;
}