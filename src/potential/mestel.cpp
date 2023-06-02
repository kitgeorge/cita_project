#include "mestel.hpp"

namespace vrs = vectors;

namespace potential {

PotentialFuncs getMestel(double v_c) {
    std::function<double(double, double, double)>
    potential = [=] (double R, double phi, double t) {
        return pow(v_c, 2)*log(R);
    };
    std::function<std::array<double, 2>(double, double, double)>
    polar_force = [=] (double R, double phi, double t) {
        double f_R = -pow(v_c, 2)/R;
        std::array<double, 2> output = {f_R, 0};
        return output;
    };
    potential::PotentialFuncs output(potential, polar_force);
    return output;
}

}