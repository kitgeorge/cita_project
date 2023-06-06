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
        double gamma = 2*pot.Omega(R_E)/pot.kappa(R_E);
        return gamma*surface_density(R_E)/(2*std::numbers::pi*pow(sigma_R(R_E), 2))
                *exp(pot.Omega(R_E)*(L - pot.LcGivenRc(R_E))/pow(sigma_R(R_E), 2));
    };
    return output;
}

}