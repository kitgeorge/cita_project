#include "potential_from_density.hpp"
#include "axsym_funcs.hpp"
#include "add_functions.hpp"
#include "integrate_rk4.hpp"
#include "nd_vectors.hpp"
#include <memory>

namespace n_body {

/**
 * Evolves coordinates of masses under the action of their BFE potential
 *
 *
 *
 *
 */
class BFENBody {
    const std::shared_ptr<const basis_functions::BFE> expansion;

    const double timestep;
    const int save_interval;
    const double N_timesteps;
    const int N_particles;
    
    // Mass of each (N_particles)
    const std::vector<double> masses;
    
    const potential::AxsymFuncs background;


    // Temporary coordinates of each particle (N_particles)
    std::vector<std::array<std::array<double, 2>, 2>>
    coords;

    // Temporary PotentialFromDensity from temporary coords
    // Optional so that it can be destroyed and a new object constructed,
    // so const-correctness within each object is ok. Perhaps we should
    // refactor so that these things aren't constant within 
    // PotentialFromDensity
    basis_functions::PotentialFromDensity bfe_pot;
    const potential::PotentialFuncs init;

    // Tabulated coordinates ((N_timesteps/save_interval + 1) * N_particles)
    utility::vector2d<std::array<std::array<double, 2>, 2>>
    saved_trajectories;

    // Tabulated BFE coefficients for all timesteps (N_timesteps + 1)
    // (N_timesteps timesteps + initial values)
    utility::vector3d<std::complex<double>>
    bfe_coefficients;

    // Tabulated norms of BFE coefficients for all timesteps
    std::vector<double> bfe_coefficient_norms;

    void getPotential();

    /**
     * This is messy but a quick fix to preexisting messiness.
     *
     * We can't construct bfe_pot with getPotential(), as the latter
     * doesn't return something but directly modifies the member.
     * Therefore, to construct the const init, which depends on 
     * bfe_pot, we have to invoke a function which, as a side
     * effect, first constructs bfe_pot.
     */
    potential::PotentialFuncs getInit();

    std::vector<std::array<std::array<double, 2>, 2>>
    iterate();

    public:
        /**
         * Constructor, which carries out simulation
         *
         * 
         *
         * @param timestep_ simulation timestep
         * @param save_interval_ trajectory is saved on one out every 
         * save_interval timesteps
         * @param integration_time Total time to simulate
         * @param background_ Background potential (eg. Mestel). This includes
         * the a contribution from the unperturbed DF
         * @param masses_ ordered list of particle masses
         * @param init_coords ordered list of intial particle polar coordinates
         * in phase space: {{R, phi}, {v_R, v_phi}}
         *
         * 
         *
         */
        BFENBody(double timestep_, int save_interval_, 
                 double integration_time, int N_particles_,
                 const potential::AxsymFuncs& background_,
                 const std::vector<double>& masses_,
                 const std::vector<std::array<std::array<double, 2>, 2>>&
                 init_coords);
        
        BFENBody(const BFENBody& old);

        /**
         * Retrieve saved particle trajectories.
         *
         * @return an (N_timesteps/save_interval + 1)*N_particles 
         * vector (indexed as (timestep, particle index)) of phase space
         * coordinates {{R, phi}, {v_R, v_phi}}
         */
        utility::vector2d<std::array<std::array<double, 2>, 2>>
        getTrajectories();

        /**
         * Retrieve BFE coefficients of density/potential for each timestep
         *
         * @return an (N_timesteps + 1)*(n_max + 1)*(l_max + 1) vector
         * of BFE coefficients, indexed as (timestep, n, l). BFE coefficient
         * for negative l is complex conjugate of positive-l counterpart*/
        utility::vector3d<std::complex<double>>
        getBFECoefficients();

        /**
         * Retrieve norm of BFE coefficients (by inner product) for each timestep
         *
         * Variations in the norm over time could mean that truncation of the 
         * BFE is significant enough to introduce errors.
         *
         * @return an N_timesteps + 1 vector, giving the norm for each timestep
         **/
        std::vector<double>
        getBFECoefficientNorms();

};

}
