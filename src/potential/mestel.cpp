#include "mestel.hpp"

namespace vectors = vrs;

namespace potential {

PotentialFuncs getMestel(double v_c) {
    std::function<double(double, double, double)>
    potential = [=] (double R, double phi, double t) {
        return pow(v_c, 2)*log(R);
    };
    std::function<vrs::Force(double, double, double)>
    force = [=] (double R, double phi, double t) {
        double f_R = -pow(v_c, 2)/R;
        vrs::Force output({f_R, 0, 0});
        return output;
    };
    potential::PotentialFuncs output(potential, force);
    return output;
}

}