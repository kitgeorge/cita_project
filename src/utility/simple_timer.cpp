#include "simple_timer.hpp"

namespace utility {

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

std::chrono::duration<double> SimpleTimer::duration() const {
    return duration.value();
}

}