#include <functional>
#include <numbers>
#include <cmath>
#include <array>
#include "axsym_funcs.hpp"
#include "mapXVtoAA2D.hpp"
#include "units.hpp"

namespace actions {

/**
 * Gets the tanh-sinh quadrature integrand for calculating T and theta_R
 *
 * Calculating the time period T of radial oscillations, and/or 
 * calculating the radial angle of an R and p_R during an oscillation,
 * involves integrating the motion in R between pericentre and apocentre.
 * The standard integrand diverges at the pericentre/apocentre, so
 * we change variables to integrate over some x between +/- 1. We then change
 * variables again, to x = tanh(pi/2 * sinh(u)), and integrate over all u. 
 * This is called tanh-sinh quadrature, and this function gives us the integrand
 * as a function of u, given the potential and parameters of the orbit
 *
 * @param E E of orbit
 * @param L L of orbit
 * @param pot Unperturbed potential
 *
 * @return integrand function, which can be usedfor calculating T and/or 
 * radial angle
 */
std::function<double(double)> 
getThetaRSTQIntegrand(potential::AxsymFuncs pot,
                      double E, double L);

/**
 * Uses Rimpei's code to find apsis for an orbit
 *
 * @param E E of orbit
 * @param L L of orbit
 * @param pot Unperturbed potential in which to find Apsis
 * givn E and L of orbit
 * 
 * @return {pericentre radius, apocentre radius}
 */
std::array<double, 2>
findApsis(potential::AxsymFuncs pot, double E, double L);

/**
 * Calculates u for a given R value of an orbit
 * 
 * @param E E of orbit
 * @param L L of orbit
 * @param pot Unperturbed potential
 * @param R R value for which to calculate u
 *
 * @return u value for R
 *
 * @note this is useful for calculating theta_R for a point between
 * the apses. That will also depend on the sign of v_R.
 */
double calculate_u(const potential::AxsymFuncs& pot, double E, double L, double R);

/**
 * Get effective potential function of unperturbed potential given L
 * @param L L of orbit
 * @param pot Unperturbed potential 
 *
 * @return Effective potential in R as an std::function
 */
std::function<double(double)>
getPhiEff(potential::AxsymFuncs pot, double L);

}