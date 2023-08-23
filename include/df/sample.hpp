#include <array>
#include <cstdlib>
#include <functional>
#include <numbers>
#include "theta_r_integrator.hpp"

namespace df {

std::array<std::array<double, 2>, 2> 
getDFSample(std::function<double(double, double)> df, 
                          double R_min, double R_max, 
                          double v_R_max, double delta_v_phi_max, double v_c);

/**
 * Randomly samples a DF in E and L
 *
 * @param df f(E, L). Must not be greater than 1 anywhere within bounds
 * @param bounds range over which to sample:
 * {{E_min, E_max}, {L_min, L_max}}.
 * We sample E_min <= E < E_max, L_min <= L < L_max;
 *
 * @return a randomly sampled pair {E, L}
 *
 * @note the closer the DF is to a normalisation where its
 * maximum value is 1, the fewer iterations it is expected to
 * take to find a sample. The function works by choosing a
 * random {E, L} within the bounds, generating a random number
 * (< 1) and accepting the {E, L} if that number is smaller than the 
 * df value (repeating otherwise).
 */
std::array<double, 2>
getDFSampleEL(std::function<double(double, double)> df,
              std::array<std::array<double, 2>, 2> bounds);

/**
 * Randomly samples a DF f(E, L) and gives phase space coordinates
 *
 * Randomly samples in E and L, then uniformly samples in angles
 * for an orbit with the actions corresponding to E and L in the 
 * given unperturbed potential.
 * Returns the corresponding phase space coordinates.
 *
 * @param df f(E, L)
 * @param E_L_bounds sampling range for (E, L): 
 * {{E_min, E_max}, {L_min, L_max}}. We sample E_min <= E < E_max,
 * L_min <= L < L_max
 * @param axsym_potential Unperturbed potential in which to find 
 * orbit for this (E, L)
 * @param u_max defines integration range in u for calculating radial 
 * coordinates from radial angle, when uniformly sampling in angles
 * @param N_u_intervals Number of integration intervals in u, see
 * ThetaRIntegrator 
 * @param N_u_iterate Number of iterations for calculating R from theta_R,
 * see ThetaRIntegrator
 *
 * @return sample {{R, phi}, {p_R, p_phi}}
 *
 * @note considering this uses df, E_L_bounds and axsym_potential, 
 * perhaps it should be absorbed into the taperedDF class? Unless
 * perhaps we want to use it for other such classes. Hmm. Maybe
 * we'll want it to act on a taperedDF object, or a parent object
 * if taperedDF ends up with siblings.
 */
std::array<std::array<double, 2>, 2>
getDFSampleViaEL(std::function<double(double, double)> df,
                 std::array<std::array<double, 2>, 2> E_L_bounds,
                 potential::AxsymFuncs axsym_potential, 
                 double u_max, int N_u_intervals, int N_u_iterate);
}