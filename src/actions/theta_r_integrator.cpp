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

double ThetaRIntegrator::calculateThetaR(std::array<double, 2> R_coords) {
    double u = calculate_u(pot, E, L, R_coords[0]);
    double interval = 2*u_max/N_intervals;
    double N_intervals_u = (int)((u - (-u_max))/interval);
    double integral = 0;
    for(int i = 0; i < N_intervals_u; ++i) {
        integral += (integrand(-u_max + i*interval)
                     + integrand(-u_max + (i + 1)*interval))
                     /2*interval;
    }
    // I'm pretty sure the integrand is doubled so that the integral only goes
    // out and not back (I'm even more sure after checking calculateR);
    double T = getT();
    double theta_R = integral/T*std::numbers::pi;
    if(R_coords[1] < 0) {
        theta_R = 2*std::numbers::pi - theta_R;
    }
    // Apses are calculated in Rimpei's code. I believe that, within some small 
    // precision, the apses may be wider than allowed by E (unless Rimpei's code
    // forces them to be narrower, which doesn't seem to be the case). Therefore
    // these integrals can become undefined at the apses, so to avoid this we don't
    // go above a certain u (so we stay a certain distance from the apses).
    // Let's try u_max = 3.
    if(std::isnan(theta_R)) {
        std::cout << "NAN theta_R: " << E << ", " << L << ", "
                  << u << ", " << integral << ", " << T << std::endl;
    }
    return theta_R;
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