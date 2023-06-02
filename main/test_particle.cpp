#include <numbers>
#include "flatten.hpp"
#include "vector_io.hpp"
#include "add_functions.hpp"
#include "units.hpp"
#include "coords.hpp"
#include "mestel.hpp"
#include "spiral.hpp"
#include "sample.hpp"
#include "integrate_rk4.hpp"
#include "execute_in_parallel.hpp"

namespace vrs = vectors;
namespace ptl = potential;

int main() {
    int N_test_particles = 36;
    int N_timesteps = 1e5;

    double mestel_vc = 220*consts::km;
    double R_0 = 8*consts::kpc;
    double spiral_amplitude_fraction = 0.0001;
    double local_pitch_angle = std::numbers::pi/6;
    int m = 2;
    double k_R = m/R_0/tan(local_pitch_angle);
    double local_Omega = mestel_vc/R_0;
    double local_kappa = local_Omega*sqrt(2);
    double pattern_speed = local_Omega + local_kappa/m;

    double integration_time = 20*2*std::numbers::pi/local_Omega;
    double timestep = integration_time/N_timesteps;

    ptl::PotentialFuncs axsym_potential = ptl::getMestel(mestel_vc);
    ptl::PotentialFuncs spiral_potential 
    = ptl::getSpiralPotential(m, k_R, spiral_amplitude_fraction
                                         *axsym_potential.potential(R_0, 0, 0),
                              pattern_speed);
    ptl::PotentialFuncs total_potential(axsym_potential, spiral_potential);
    // ptl::PotentialFuncs total_potential = axsym_potential;
    

    double gamma = 2*local_Omega/local_kappa;
 

    std::vector<vrs::Coords2d>
    initial_conditions(N_test_particles);
    for(int i = 0; i < N_test_particles; ++i) {
        initial_conditions[i].setPolar({{ {{R_0, (double)i/18*std::numbers::pi}},
                                          {{30*consts::km, mestel_vc}} }});
    }
    






    std::vector<std::function<std::vector<vrs::Coords2d>()>>
    tp_integration_functions(N_test_particles); 
    for(int i = 0; i < N_test_particles; ++i) {
        tp_integration_functions[i] 
        = tp_integration::getTpIntegrationFunction(total_potential, initial_conditions[i], 
                                                   0, timestep, N_timesteps);
    }

    std::vector<std::vector<vrs::Coords2d>>
    trajectories = multithreading::executeInParallel(tp_integration_functions);
    
    // std::vector<std::vector<vrs::Coords2d>>
    // sampled_trajectories(1000);
    // for(int i = 0; i < 1000; ++i) {
    //     sampled_trajectories[i] = trajectories[i*1000];
    // }
    

    utility::writeCsv("../data/test_particle_1.csv", vrs::getPolarsFlat(utility::flatten(trajectories)));

}

    