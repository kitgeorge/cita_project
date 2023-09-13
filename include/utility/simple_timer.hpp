#include <chrono>
#include <cassert>
#include <optional>

namespace utility {

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

}