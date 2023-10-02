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

std::array<std::array<double, 2>, 2>
getDFSampleViaEL(std::function<double(double, double)> df,
                 std::array<std::array<double, 2>, 2> E_L_bounds,
                 potential::AxsymFuncs axsym_potential, 
                 double u_max, int N_u_intervals, int N_u_iterate) {
    std::array<double, 2>
    sample_E_L = df::getDFSampleEL(df, E_L_bounds);
    double angle_normalisation = (double)RAND_MAX/(2*std::numbers::pi);
    double theta_R = (double)std::rand()/angle_normalisation;
    std::cout << "theta_R: " << theta_R << std::endl;
    actions::ThetaRIntegrator integrator(axsym_potential, sample_E_L[0],
                                         sample_E_L[1], u_max,
                                         N_u_intervals, N_u_iterate);
    std::array<std::array<double, 2>, 2> output;
    std::array<double, 2> R_coords = integrator.getCoords(theta_R);
    output[0][0] = R_coords[0];
    output[1][0] = R_coords[1];
    output[0][1] = (double)std::rand()/angle_normalisation;
    output[1][1] = sample_E_L[1]/output[0][0];
    return output;    
}


}