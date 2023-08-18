#include "tapered_df.hpp"

namespace df {

TaperedDF::TaperedDF(double v_c_, double R_0_, double active_fraction,
                     std::array<double, 2> taper_radii_, 
                     std::array<double, 2> taper_indices_,
                     std::array<double, 2> cutoff_radii,
                     const std::function<double(double)>& target_Q):
    v_c(v_c_), R_0(R_0_), mestel_potential(potential::getMestel(v_c, R_0)),
    E_bounds(getEBounds(cutoff_radii)), L_bounds(getLBounds(cutoff_radii)),
    dehnen_df(calculateDehnenDF(active_fraction, target_Q)),
    taper_radii(taper_radii_), taper_indices(taper_indices_),
    tapered_df(calculateTaperedDF()) {}

TaperedDF::TaperedDF(const TaperedDF& old):v_c(old.v_c), R_0(old.R_0),
    mestel_potential(old.mestel_potential),
    E_bounds(old.E_bounds), L_bounds(old.L_bounds), dehnen_df(old.dehnen_df),
    taper_radii(old.taper_radii), taper_indices(old.taper_indices),
    tapered_df(old.tapered_df) {}

std::array<double, 2>
TaperedDF::getEBounds(std::array<double, 2> cutoff_radii) const {
    std::array<double, 2> output;
    for(int i = 0; i < 2; ++i) {
        output[i] = mestel_potential.EcGivenRc(cutoff_radii[i]);
    }
    return output;
}

std::array<double, 2>
TaperedDF::getLBounds(std::array<double, 2> cutoff_radii) const {
    std::array<double, 2> output;
    for(int i = 0; i < 2; ++i) {
        output[i] = mestel_potential.LcGivenRc(cutoff_radii[i]);
    }
    return output;
}

std::function<double(double, double)>
TaperedDF::calculateDehnenDF(double active_fraction,
                             const std::function<double(double)>& target_Q)
                             const {
    std::function<double(double)> 
    target_density = [this, active_fraction] (double R) {
        return active_fraction*pow(v_c, 2)/(2*std::numbers::pi*Units::G*R);
    };
    std::function<double(double)>
    target_sigma = [this, target_density, target_Q] (double R) {
        return target_Q(R)*3.36*Units::G*target_density(R)
               /mestel_potential.kappa(R);
    };
    return getNormDehnenDF(mestel_potential, target_density,
                           target_sigma, E_bounds);
}

double TaperedDF::innerTaper(double L) const {
    double numerator = pow(L, taper_indices[0]);
    double denominator = pow(taper_radii[0]*v_c, taper_indices[0])
                         + numerator;
    return numerator/denominator;
}

double TaperedDF::outerTaper(double L) const {
    double denominator = 1 + pow(L/(taper_radii[1]*v_c), taper_indices[1]);
    return 1/denominator;
}

double tdfFunction(double E, double L) const {
    return innerTaper(L)*outerTaper(L)*dehnen_df(E, L);
}

std::function<double(double, double)>
TaperedDF::getTaperedDF() const {
    return [*this] (double E, double L) {
        return tdfFunction(E, L);
    };
} 

}