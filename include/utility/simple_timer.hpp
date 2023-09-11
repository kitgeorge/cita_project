#include <chrono>
#include <cassert>

namespace utility {

class SimpleTimer {
    std::optional<const std::chrono::time_point<std::chrono::steady_clock>> 
    start_time;

    std::optional<const std::chrono::duration<double>> duration;

    public:
        void start();
        void stop();
        std::chrono::duration<double> getDuration() const;
};

}