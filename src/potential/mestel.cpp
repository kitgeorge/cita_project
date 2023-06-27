#include "mestel.hpp"


namespace vrs = vectors;

namespace potential {

AxsymFuncs getMestel(double v_c, double R_0) {
    std::function<double(double, double, double)>
    potential = [=] (double R, double phi, double t) {
        return pow(v_c, 2)*log(R/R_0);
    };
    std::function<std::array<double, 2>(double, double, double)>
    polar_force = [=] (double R, double phi, double t) {
        double f_R = -pow(v_c, 2)/R;
        std::array<double, 2> output = {f_R, 0};
        return output;
    };
    std::function<double(double)>
    dPhidR = [v_c] (double R) {
        return pow(v_c, 2)/R;
    };
    std::function<double(double)>
    d2PhidR2 = [v_c] (double R) {
        return -pow(v_c/R, 2);
    };
    std::function<double(double)>
    RcGivenL = [v_c] (double L) {
        return L/v_c;
    };
    std::function<double(double)>
    RcGivenE = [v_c, R_0] (double E) {
        return R_0*exp(E/pow(v_c, 2) - 0.5);
    };
    std::function<double(double)>
    RcGivenOmega = [v_c] (double Omega) {
        return v_c/Omega;
    };
    std::function<double(double)>
    vcGivenRc = [v_c] (double R_c) {
        return v_c;
    };
    std::function<double(double)>
    LcGivenRc = [v_c] (double R_c) {
        return v_c*R_c;
    };
    std::function<double(double)>
    EcGivenRc = [v_c, potential] (double R_c) {
        return pow(v_c, 2)/2 + potential(R_c, 0, 0);
    };
    std::function<double(double)>
    EcGivenL = [v_c, potential] (double L) {
        return pow(v_c, 2)/2 + potential(L/v_c, 0, 0);
    };
    std::function<double(double)>
    Omega = [v_c] (double R) {
        return v_c/R;
    };
    std::function<double(double)>
    kappa = [Omega] (double R) {
        return sqrt(2)*Omega(R);
    };
    std::function<double(int, int, double)>
    resonanceRadius = [v_c] (int N_R, int N_phi, double Omega_p) {
        return (1 + sqrt(2)*(double)N_R/N_phi)*v_c/Omega_p;
    };
    AxsymFuncs output(potential, polar_force, 
    {{dPhidR, d2PhidR2, RcGivenL, RcGivenE, RcGivenOmega, vcGivenRc, LcGivenRc,
      EcGivenRc, EcGivenL, Omega, kappa}}, resonanceRadius);
    return output;
}

}