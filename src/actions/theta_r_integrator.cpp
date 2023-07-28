#include "theta_r_integrator.hpp"

namespace actions {



ThetaRIntegrator::ThetaRIntegrator(potential::AxsymFuncs pot_, 
                                   double E_, double L_,
                                   double u_max_, int N_intervals_,
                                   int N_iterate_): 
                                   pot(pot_),
                                   apsis(findApsis(pot_, E_, L_)),
                                   u_max(u_max_), N_intervals(N_intervals_),
                                   N_iterate(N_iterate_), E(E_), L(L_), T(0) {
    integrand = getThetaRSTQIntegrand(pot, E, L);
    calculateT();
}

void ThetaRIntegrator::calculateT() {
    // Use Trapezium rule for integration in u (symmetry)
    double interval = 2*u_max/N_intervals;
    for(int i = 0; i < N_intervals; ++i) {
        T += (integrand(-u_max + i*interval) 
              + integrand(-u_max + (i + 1)*interval))
             /2*interval;
    }
}

double ThetaRIntegrator::getT() {
    return T;
}

double ThetaRIntegrator::getRFromu(double u) {
    double x = tanh(std::numbers::pi/2*sinh(u));
    return (apsis[0] + apsis[1])/2 + x*(apsis[1] - apsis[0])/2;
}

double ThetaRIntegrator::calculateR(double theta_R) {
    theta_R = fmod(theta_R, 2*std::numbers::pi);
    if(theta_R > std::numbers::pi) {
        theta_R = 2*std::numbers::pi - theta_R; 
    }
    std::array<double, 2> bounds = {-u_max, u_max};
    double theta_integral = 0;


    double sub_interval;

    for(int i = 0; i < N_iterate; ++i) {
        sub_interval = (bounds[1] - bounds[0])/N_intervals;
        for(int j = 0; j < N_intervals; ++j) {
            double increment = (integrand(bounds[0] + j*sub_interval)
                                + integrand(bounds[0] + (j + 1)*sub_interval))
                                /2*sub_interval*std::numbers::pi/T;
            if(theta_integral + increment > theta_R) {
                bounds = {bounds[0] + j*sub_interval,
                            bounds[0] + (j + 1)*sub_interval};
                break;
            }
            else {
                theta_integral += increment;
            }

        }
    }
    double u = bounds[0];
    return getRFromu(u); 
}

std::array<double, 2>
ThetaRIntegrator::getCoords(double theta_R) {
    std::function<double(double)>
    Phi_eff = getPhiEff(pot, L);
    double R = calculateR(theta_R);
    double p_R = sqrt(2*(E - Phi_eff(R)))*(2*(theta_R <= std::numbers::pi) - 1);
    std::array<double, 2> output = {R, p_R};
    return output;
}

}