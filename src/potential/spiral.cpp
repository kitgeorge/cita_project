#include "spiral.hpp"

namespace potential {

PotentialFuncs getSpiralPotential(double m, double k_R, 
                                  double amplitude,
                                  double pattern_speed) {

    std::function<double(double, double, double)>
    potential = [=] (double R, double phi, double t) {
        return amplitude*cos(m*(phi - pattern_speed*t) + k_R*R);
    };

    std::function<vrs::Force(double, double, double)>
    force = [=] (double R, double phi, double t) {
        double f_R = k_R*amplitude
                     *sin(m*(phi - pattern_speed*t) + k_R*R);
        double f_phi = m*amplitude
                     *sin(m*(phi - pattern_speed*t) + k_R*R);
        vrs::Force output({f_R, f_phi, 0});
        return output;
    };
    
    potential::PotentialFuncs output(potential, force);
    return output;
}

}