#include "theta_r_integrand.hpp"
#include <iostream>


namespace actions {

class ThetaRIntegrator {
    potential::AxsymFuncs pot;
    std::function<double(double)> integrand;
    std::array<double, 2> apsis;
    double u_max;
    int N_intervals;
    int N_iterate;
    double T;
    double E;
    double L;

    public:
        ThetaRIntegrator(potential::AxsymFuncs pot_, 
                         double E_, double L_,
                         double u_max_, int N_intervals_, 
                         int N_iterate_);
    void calculateT();
    double getT();
    double getRFromu(double u);
    double calculateR(double theta_R);
    std::array<double, 2> getCoords(double theta_R);

};

}