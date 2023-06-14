#include "spiral.hpp"

namespace potential {

PotentialFuncs getSpiralPotential(double m, double k_R, 
                                  double amplitude,
                                  double pattern_speed,
                                  double initial_phase) {

    std::function<double(double, double, double)>
    potential = [=] (double R, double phi, double t) {
        return amplitude*cos(m*(phi - pattern_speed*t) + k_R*R
                             + initial_phase);
    };

    std::function<std::array<double, 2>(double, double, double)>
    polar_force = [=] (double R, double phi, double t) {
        double f_R = k_R*amplitude
                     *sin(m*(phi - pattern_speed*t) + k_R*R
                          + initial_phase);
        double f_phi = m/R*amplitude
                     *sin(m*(phi - pattern_speed*t) + k_R*R
                          + initial_phase);
        std::array<double, 2> output = {f_R, f_phi};
        return output;
    };
    
    potential::PotentialFuncs output(potential, polar_force);
    return output;
}

// PotentialFuncs getLogSpiralPotential(int m, double pitch_angle,
//                                      double amplitude,
//                                      double pattern_speed,
//                                      double R_0, double initial_phase) {
//     std::function<double(double, double, double)>
//     potential = [=] (double R, double phi, double t) {
//         return amplitude*cos(m*(phi - pattern_speed*t)
//                              + m/tan(pitch_angle)*ln(R/R_0)
//                              + initial_phase);
//     };
//     std::function<std::array<double, 2>(double, double, double)>
//     polar_force = [=] (double R, double phi, double t) {
//         double wave = amplitude*sin(m*(phi - pattern_speed*t)
//                              + m/tan(pitch_angle)*ln(R/R_0)
//                              + initial_phase);
//         double k_phi = m/R;
//         double k_R = k_phi/tan(pitch_angle);
//         std::array<double, 2> output = {k_R*wave, k_phi*wave};
//         return output;
//     }
//     potential::PotentialFuncs output(potential, polar_force);
//     return output;
// }

// PotentialFuncs getSpiralPacket(AxsymFuncs pot, double R_0,
//                                )


}