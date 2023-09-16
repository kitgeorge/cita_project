#include <cstdlib>
#include <optional>
#include <cassert>
#include <vector>
#include <iostream>
#include <chrono>

class SimpleTimer {
    std::optional<const std::chrono::time_point<std::chrono::steady_clock>> 
    start_time;

    std::optional<const std::chrono::duration<double>> duration;

    public:
        void start();
        void stop();
        std::chrono::duration<double> getDuration() const;
        double getDuration_ms() const;
        double getDuration_us() const;
        double getDuration_ns() const;
};


void SimpleTimer::start() {
    start_time.reset();
    duration.reset();
    start_time.emplace(std::chrono::steady_clock::now());
}

void SimpleTimer::stop() {
    assert(start_time);
    duration.emplace(std::chrono::steady_clock::now() - start_time.value());
    start_time.reset();
}

std::chrono::duration<double> SimpleTimer::getDuration() const {
    return duration.value();
}

double SimpleTimer::getDuration_ms() const {
    return std::chrono::duration_cast<std::chrono::milliseconds>(duration.value()).count();
}

double SimpleTimer::getDuration_us() const {
    return std::chrono::duration_cast<std::chrono::microseconds>(duration.value()).count();
}

double SimpleTimer::getDuration_ns() const {
    return std::chrono::duration_cast<std::chrono::nanoseconds>(duration.value()).count();
}


int main() {
    SimpleTimer timer;
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
        // std::cout << b[i];
    }
    std::cout << std::endl << "Read double: " << timer.getDuration_ns()/N_a
              << "ns" << std::endl;
}