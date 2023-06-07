#include "dehnen_df.hpp"

namespace df {

std::function<double(double, double)>
getDehnenDF(const potential::AxsymFuncs& pot, 
            std::function<double(double)> surface_density,
            std::function<double(double)> sigma_R) {
    // We assume that R_E (=R_c(E)) increases with E
    std::function<double(double, double)>
    output = [=] (double E, double L) {
        double R_E = pot.RcGivenE(E);
        // Eliminate invalid (E, L) combination (so we can straightforwardly 
        // tabulate in (E, L))
        if(L > pot.LcGivenRc(R_E)) {
            return 0.0;
        }
        double gamma = 2*pot.Omega(R_E)/pot.kappa(R_E);
        return gamma*surface_density(R_E)/(2*std::numbers::pi*pow(sigma_R(R_E), 2))
                *exp(pot.Omega(R_E)*(L - pot.LcGivenRc(R_E))/pow(sigma_R(R_E), 2));
    };
    return output;
}

std::function<double(double, double)>
getNormDehnenDF(const potential::AxsymFuncs& pot,
                std::function<double(double)> surface_density,
                std::function<double(double)> sigma_R,
                std::array<double, 2> E_bounds) {
    // Normalise so that f_max = 1 within (E, L) bounds
    std::function<double(double)>
    prefactor = [=] (double E) {
        double R_E = pot.RcGivenE(E);
        double gamma = 2*pot.Omega(R_E)/pot.kappa(R_E);
        return gamma*surface_density(R_E)/(2*std::numbers::pi*pow(sigma_R(R_E), 2));
    };
    double norm = std::max(prefactor(E_bounds[0]), prefactor(E_bounds[1]));
    std::cout << norm << std::endl;
    return utility::multiplyFunction(getDehnenDF(pot, surface_density, sigma_R),
                                     1/norm); 
}

}