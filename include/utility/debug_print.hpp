#include <iostream>
#include <string>
#include <array>
#include <mutex>
#pragma once

namespace utility {

static const int N_debug_channels = 10;
extern std::array<bool, N_debug_channels> channel_open;

void open_channel(int channel);
void close_channel(int channel);

void debug_print(const std::string& statement, int channel);

}