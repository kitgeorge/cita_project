#include "debug_print.hpp"

namespace utility {

// Set all debug channels closed by default
std::array<bool, N_debug_channels> 
channel_open = [N_debug_channels]() {
    std::array<bool, N_debug_channels> output;
    for(int i = 0; i < N_debug_channels; ++i) {
        output[i] = 0;
    }
    return output;
}();

void open_channel(int channel) {
    channel_open[channel] = 1;
}

void close_channel(int channel) {
    channel_open[channel] = 0;
}

void debug_print(const std::string& statement, int channel) {
    if(channel_open[channel]) {
        std::cout << statement << std::endl;
    }
}




}