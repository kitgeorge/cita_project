#include "potential_funcs.hpp"

namespace potential {

void PotentialFuncs::multiply(std::function<double(double, double, double)>
                              envelope) {
    std::function<double(double, double, double)>
    potential_0 = potential;
    std::function<std::array<double, 2>(double, double, double)>
    polar_force_0 = polar_force;
    std::function<std::array<double, 2>(double, double, double)>
    cartesian_force_0 = cartesian_force;

    potential = [potential_0, envelope] (double R, double phi, double t) {
        return envelope(R, phi, t)*potential_0(R, phi, t);
    };

    polar_force = [envelope, polar_force_0] (double R, double phi, double t) {
        std::array<double, 2> f = polar_force_0(R, phi, t);
        for(int i = 0; i < 2; ++i) {
            f[i] = envelope(R, phi, t)*f[i];
        }
        return f;
    };

    cartesian_force = [envelope, cartesian_force_0](double x, double y, double t) {
        std::array<double, 2> f = cartesian_force_0(x, y, t);
        double R = sqrt(pow(x, 2) + pow(y, 2));
        double phi = atan(y/x) + std::numbers::pi*(x < 0);
        for(int i = 0; i < 2; ++i) {
            f[i] = envelope(R, phi, t)*f[i];
        }
        return f;
    };

    }
} 