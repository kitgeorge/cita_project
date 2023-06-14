#include "sample.hpp"
#include <iostream>

namespace df {

std::array<std::array<double, 2>, 2> 
getDFSample(std::function<double(double, double)> df, 
                          double R_min, double R_max, 
                          double v_R_max, double delta_v_phi_max,
                          double v_c) {
    // df must have maximum value <= 1
    double R, phi, v_R, v_phi;
    int accepted_flag = 0;
    while(accepted_flag == 0) {
        R = R_min + (double)std::rand()/RAND_MAX*(R_max - R_min);
        phi = (double)std::rand()/RAND_MAX*2*std::numbers::pi;
        v_R = (2*(double)std::rand()/RAND_MAX - 1)*v_R_max;
        v_phi = v_c + (2*(double)std::rand()/RAND_MAX - 1)*delta_v_phi_max;
        if(df(v_R, v_phi) > (double)std::rand()/RAND_MAX) {
            accepted_flag = 1;
        }
    }
    std::array<std::array<double, 2>, 2>
    output = {{ {{R, phi}},
                {{v_R, v_phi}} }};
    return output;
}

std::array<double, 2>
getDFSampleEL(std::function<double(double, double)> df,
              std::array<std::array<double, 2>, 2> bounds) {
    // Should normalise df such that its maximal value over sampled bounds is 1
    // (any lower and this function will be slower, any higher and it saturates)
    double E, L;
    int accepted_flag = 0;
    while(accepted_flag == 0) {
        E = bounds[0][0] + (double)std::rand()/RAND_MAX
                                              *(bounds[0][1] - bounds[0][0]);
        L = bounds[1][0] + (double)std::rand()/RAND_MAX
                                              *(bounds[1][1] - bounds[1][0]);
        // std::cout << E << ", " << L << ", " << df(E, L) << std::endl;
        if(df(E, L) > (double)std::rand()/RAND_MAX) {
            accepted_flag = 1;
        }
    }
    std::array<double, 2> output = {E, L};
    return output;
}


}