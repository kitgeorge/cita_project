#include "theta_r_integrand.hpp"

namespace actions {

std::function<double(double)> 
getThetaRSTQIntegrand(potential::AxsymFuncs pot,
                      double E, double L) {
    std::function<double(double)>
    Phi_eff = getPhiEff(pot, L);
    
    std::array<double, 2> 
    apsis = findApsis(pot, E, L);
    double R_apo = apsis[1];
    double R_peri = apsis[0];
    std::function <double(double)>
    R_u = [=] (double u) {
        return (R_apo + R_peri)/2 
                + (R_apo - R_peri)/2*tanh(std::numbers::pi/2*sinh(u));
    };
    std::function <double(double)>
    jacobian = [=] (double u) {
        return (R_apo - R_peri)/2*pow(1/cosh(std::numbers::pi/2*sinh(u)), 2)
               *std::numbers::pi/2*cosh(u);
    };
    return [=] (double u) {
        return sqrt(2)*jacobian(u)/sqrt(E - Phi_eff(R_u(u)));
    };

}

double calculate_u(const potential::AxsymFuncs& pot, double E, double L, double R) {
    std::array<double, 2> apsis = findApsis(pot, E, L);
    double R_0 = (apsis[0] + apsis[1])/2;
    double x = (R - R_0)/(0.5*(apsis[1] - apsis[0]));
    double u = std::asinh(2/std::numbers::pi*std::atanh(x));
    return u;
}


std::function<double(double)>
getPhiEff(potential::AxsymFuncs pot, double L) {
    return [=] (double R) {
        return pow(L, 2)/(2*pow(R, 2)) + pot.potential_R(R);
    };
}


std::array<double, 2>
findApsis(potential::AxsymFuncs pot, double E, double L) {
    RC::MapXVtoAA2D finding_apsis(&pot);
    finding_apsis.findApsis(L, E);
    std::array<double, 2> output;
    output = {finding_apsis.r_peri, finding_apsis.r_apo};
    return output;
}

}