#include "sample.hpp"

namespace df {

std::array<std::array<double, 2>, 2> 
getDFSample(std::function<double(double, double)> df, 
                          double R_min, double R_max, 
                          double v_R_max, double v_phi_max) {
    // df must have maximum value <= 1
    double R, phi, v_R, v_phi;
    int accepted_flag = 0;
    while(accepted_flag == 0) {
        R = R_min + (double)std::rand()/RAND_MAX*(R_max - R_min);
        phi = (double)std::rand()/RAND_MAX*2*std::numbers::pi;
        v_R = (2*(double)std::rand()/RAND_MAX - 1)*v_R_max;
        v_phi = (2*(double)std::rand()/RAND_MAX - 1)*v_phi_max;
        if(df(v_R, v_phi) > (double)std::rand()/RAND_MAX) {
            accepted_flag = 1;
        }
    }
    std::array<std::array<double, 2>, 2>
    output = {{ {{R, phi}},
                {{v_R, v_phi}} }};
    return output;
}

}