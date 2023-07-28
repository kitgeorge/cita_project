#include "libration_calculator.hpp"
#include <iostream>

namespace actions {

double LibrationCalculator::omega_0() const {
    return sqrt(F*G);
}

void LibrationCalculator::setICs(double J_s, double theta_s) {
    k.emplace(calculatek(J_s, theta_s));
    assert(k.value() <= 1);
    J_a = 2*sqrt(F/G)*k.value();
    C.emplace(jfuncs.cnInv((J_s - J_s_res)/J_a, pow(*k, 2)));
    if(fmod(theta_s + g + 2*std::numbers::pi, 2*std::numbers::pi) 
            > std::numbers::pi) {
        C.emplace(-1*(*C));
    }

    std::array<double, 2> check = trajectory(0);
    assert(std::abs(check[0]/theta_s - 1) < 5e-2);
    assert(std::abs(check[1]/J_s - 1) < 1e-2);
}

std::array<double, 2> LibrationCalculator::trajectory(double t) const {
    assert(k != std::nullopt);
    assert(C != std::nullopt);
    double theta_s = 2*asin(*k*jfuncs.sn(sqrt(F*G)*t + *C, pow(*k, 2))) - g;
    double J_s = J_s_res + J_a*jfuncs.cn(sqrt(F*G)*t + *C, pow(*k, 2));
    std::array<double, 2> output = {theta_s, J_s};
    return output;
}


LibrationCalculator::LibrationCalculator(double J_f_, const JresFinder& finder,
                                         const potential::LibrationFuncs& libfuncs,
                                         double delta_J_s, int N_jacobi_intervals):
                    J_f(J_f_), J_s_res(finder.calculateJsres(J_f_)),
                    jfuncs(N_jacobi_intervals), N_R(finder.N_R), N_phi(finder.N_phi),
                    F(libfuncs.F(J_f, J_s_res, N_R)), g(libfuncs.g(J_f, J_s_res, N_R)),
                    G((finder.getOmega_s(J_f, J_s_res + delta_J_s)
                      - finder.getOmega_s(J_f, J_s_res - delta_J_s))/(2*delta_J_s)) {}                    

LibrationCalculator::LibrationCalculator(const LibrationCalculator& old):
    N_R(old.N_R), N_phi(old.N_phi), J_f(old.J_f), J_s_res(old.J_s_res),
    F(old.F), G(old.G), jfuncs(old.jfuncs), g(old.g) {}


// LibrationCalculator::LibrationCalculator(double J_f_, const JresFinder& finder,
//                                          std::function<double(double, double)> F_,
//                                          std::function<double(double, double)> g_,
//                                          double delta_J_s, int N_jacobi_intervals):
//                     J_f(J_f_), J_s_res(finder.calculateJsres(J_f_)),
//                     jfuncs(N_jacobi_intervals) {
//     N_R = finder.N_R;
//     N_phi = finder.N_phi;
//     F = F_(J_f, J_s_res);
//     g = g_(J_f, J_s_res);
//     G = (finder.getOmega_s(J_f, J_s_res + delta_J_s)
//        - finder.getOmega_s(J_f, J_s_res - delta_J_s))/(2*delta_J_s);
// }

std::array<double, 2> LibrationCalculator::separatrices(double theta_s) const {
    double J_a_0 = 2*sqrt(F/G);
    double diff = J_a_0*sqrt(1 - pow(sin((theta_s + g)/2), 2));
    std::array<double, 2> output = {J_s_res - diff,
                                    J_s_res + diff};
    return output;
}

double LibrationCalculator::calculateEp(double J_s, double theta_s) const {
    double dthetadt = G*(J_s - J_s_res);
    return pow(dthetadt, 2)/2 - F*G*cos(theta_s + g);
}

double LibrationCalculator::calculatek(double J_s, double theta_s) const {
    double E_p = calculateEp(J_s, theta_s);
    return sqrt((1 + E_p/(F*G))/2);
}

double JacobiSpecialFunctions::phi(double u, double m) const {
    // This is all a bit of a mess, should probably rewrite at some point
    // and test properly
    int sgn = 1;
    if(u < 0) {
        sgn = -1;
    }
    u = std::abs(u);
    std::function<double(double)> 
    integrand = [m] (double theta) {return 1/sqrt(1 - m*pow(sin(theta), 2));};
    // I know the function we're integrating is symmetric about pi/2, but
    // I can't be bothered dealing with the extra complexity in my code
    // and would rather take double the calculation time for this function
    // (can always come back to it)
    double dtheta = std::numbers::pi/N_intervals;
    std::vector<double> cumulative_values(N_intervals);
    double integral = 0;
    for(int i = 0; i < N_intervals; ++i) {
        cumulative_values[i] = integral;
        integral += dtheta/2*(integrand(i*dtheta) + integrand((i + 1)*dtheta));
    }
    int N_domains = static_cast<int>(u/integral);
    double remainder = std::fmod(u, integral);
    double output = N_domains*std::numbers::pi;
    for(int i = 0; i < N_intervals - 1; ++i) {
        if(remainder < cumulative_values[i + 1]) {
            double bin_fraction = (remainder - cumulative_values[i])
                                    /(cumulative_values[i + 1] - cumulative_values[i]);
            output += (i + bin_fraction)*dtheta;
            return sgn*output;
        }
    }
    double bin_fraction = (remainder - cumulative_values[N_intervals - 1])
                            /(integral - cumulative_values[N_intervals - 1]);
    output += (N_intervals - 1 + bin_fraction)*dtheta;
    return sgn*output;
}

double JacobiSpecialFunctions::cn(double u, double m) const {
    return cos(phi(u, m));
}

double JacobiSpecialFunctions::sn(double u, double m) const {
    return sin(phi(u, m));
}

double JacobiSpecialFunctions::cnInv(double x, double m) const {
    // Not sure I'm happy about the difference in precision between this and cn
    // (in how the last sub-interval of the domain is discarded here)
    std::function<double(double)>
    integrand = [m] (double theta) {return 1/sqrt(1 - m*pow(sin(theta), 2));};
    double dtheta = std::numbers::pi/N_intervals;
    if(x < -1) {
        x = -1;
    }
    if(x > 1) {
        x = 1;
    }
    int N_total = static_cast<int>(acos(x)/dtheta);
    double integral = 0;
    for(int i = 0; i < N_total; ++i) {
        integral += dtheta/2*(integrand(i*dtheta) + integrand((i + 1)*dtheta));
    }
    return integral;
}

double JacobiSpecialFunctions::K(double m) const {
    // These are all a bit dodgy, I may well come back to them if not 
    // accurate enough
    std::function<double(double)> 
    integrand = [m] (double theta) {return 1/sqrt(1 - m*pow(sin(theta), 2));};
    double dtheta = std::numbers::pi/N_intervals;
    double integral = 0;
    for(int i = 0; i < N_intervals/2; ++i) {
        integral += dtheta/2*(integrand(i*dtheta) + integrand((i + 1)*dtheta));
    }
    return integral;
}

}