#include "theta_r_integrand.hpp"
#include <iostream>


namespace actions {

/**
 * Class for calculating radial period T, and {R, p_R} given theta_R, for an orbit
 *
 * For an orbit of given E and L in an (unperturbed) potential, this class
 * calculates the period T of radial oscillations and the {R, p_R}
 * coordintates corresponding to radial angle theta_R
 *
 * @note theta_R is set to be zero at pericentre
 */
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
        /**
         * Constructor
         * 
         * @param pot_ unperturbed potential
         * @param E_ E of orbit
         * @param L_ L of orbit
         * @param u_max_ we integrate over all u and integrand which
         * tends to zero as u -> +/- infinity. We have truncate 
         * the integration range to -u_max <= u < u_max
         * @param N_intervals_ this serves two purposes. It is the number
         * of intervals in u we integrate over to find T. It is also
         * the number of intervals in u we search for the value of u 
         * (and so R) for a given theta_R. The latter is done
         * iteratively.
         * @param N_iterate_ Number of iterations for finding R
         * given theta_R
         */
        ThetaRIntegrator(potential::AxsymFuncs pot_, 
                         double E_, double L_,
                         double u_max_, int N_intervals_, 
                         int N_iterate_);
    /// Calculates and stores period T of radial oscillations for this orbit
    void calculateT();
    /// Retrieves stored T
    double getT();
    /// Calculates R from tanh-sinh quadrature integration variable u
    double getRFromu(double u);
    /**
     * Calculates R from theta_R for this orbit
     *
     * Integrating the tanh-sinh quadrature integrand over all u gives
     * T (divided by a factor 2 I believe). Integrating only up to the u value
     * corresponding to angle theta_R will give theta_R/(2 pi) * T 
     * (ignoring for now the complication of having to integrate from 
     * pericentre to apocentre and back again, which is dealt with in the
     * code). We want to find the value of R corresponding to theta_R.
     * To do this we find the corresponding value of u, by numerically integrating
     * in u only until the integrated value is theta_R/(2 pi) * T, and then 
     * taking the u value where this is reached. We then get the R value.
     * (please ignore the exact maths here as it is more complicated when 
     * you consider the two integrals, out and back. Have to dig into my notes to
     * see that).
     *
     * @param theta_R theta_R value for which to find corresponding R value
     *
     * @return corresponding R value
     *
     * @note this process is done through multiple refining iterations.
     * We start with N_intervals intervals covering all u (up to u_max). Then
     * when we find the correct interval in u, we split that interval into
     * N_intervals sub intervals, and search for the correct sub-interval of u.
     * We repeat this N_iterate times
     */
    double calculateR(double theta_R);
    /// Calculates {R, p_R} at given theta_R for this orbit
    std::array<double, 2> getCoords(double theta_R);

};

}