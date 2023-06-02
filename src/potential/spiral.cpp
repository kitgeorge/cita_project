#include "spiral.hpp"

namespace potential {

PotentialFuncs getSpiralPotential(double m, double k_R, 
                                  double amplitude,
                                  double pattern_speed) {

    std::function<double(double, double, double)>
    potential = [=] (double R, double phi, double t) {
        return amplitude*cos(m*(phi - pattern_speed*t) + k_R*R);
    };

    std::function<std::array<double, 2>(double, double, double)>
    polar_force = [=] (double R, double phi, double t) {
        double f_R = k_R*amplitude
                     *sin(m*(phi - pattern_speed*t) + k_R*R);
        double f_phi = m/R*amplitude
                     *sin(m*(phi - pattern_speed*t) + k_R*R);
        std::array<double, 2> output = {f_R, f_phi};
        return output;
    };
    
    potential::PotentialFuncs output(potential, polar_force);
    return output;
}

}