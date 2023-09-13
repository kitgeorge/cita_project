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

std::chrono::duration<double> SimpleTimer::getDuration() const {
    return duration.value();
}

double SimpleTimer::getDuration_ms() const {
    return std::chrono::duration_cast<std::chrono::milliseconds>(duration.value()).count();
}

double getDuration_us() const {
    return std::chrono::duration_cast<std::chrono::microseconds>(duration.value()).count();
}

double getDuration_ns() const {
    return std::crhono::duration_cast<std::chrono::nanoseconds>(duration.value()).count();
}

}