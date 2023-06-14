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
    int N_test_particles = 1000;
    int N_timesteps = 1e4;

    double mestel_vc = 220*Units::kms;
    double R_0 = 8;
    double spiral_amplitude_fraction = 0.1;
    double local_pitch_angle = std::numbers::pi/6;
    int m = 2;
    double k_R = m/R_0/tan(local_pitch_angle);
    double local_Omega = mestel_vc/R_0;
    double local_kappa = local_Omega*sqrt(2);
    double pattern_speed = local_Omega + local_kappa/m;

    double integration_time = 10*2*std::numbers::pi/local_Omega;
    double timestep = integration_time/N_timesteps;

    ptl::PotentialFuncs axsym_potential = ptl::getMestel(mestel_vc, R_0);
    ptl::PotentialFuncs spiral_potential 
    = ptl::getSpiralPotential(m, k_R, spiral_amplitude_fraction,
                              pattern_speed. 0);
    ptl::PotentialFuncs total_potential(axsym_potential, spiral_potential);
    

    double sigma_R = 30*consts::km;
    double gamma = 2*local_Omega/local_kappa;
    std::function<double(double, double)>
    schwarzchild_df = [=] (double v_R, double v_phi) {
        return exp(-(pow(v_R, 2) + pow(gamma*(v_phi - mestel_vc), 2))/(2*pow(sigma_R, 2)));
    };
    double R_min = 0.9*R_0;
    double R_max = 1.1*R_0;
    double v_R_max = 100*Units::kms;
    double delta_v_phi_max = v_R_max/gamma;








    std::vector<std::function<std::vector<vrs::Coords2d>()>>
    tp_integration_functions(N_test_particles); 
    for(int i = 0; i < N_test_particles; ++i) {
        vrs::Coords2d
        initial_conditions(df::getDFSample(schwarzchild_df, 
                                           R_min, R_max, 
                                           v_R_max, delta_v_phi_max, mestel_vc), 1);
        tp_integration_functions[i] 
        = tp_integration::getTpIntegrationFunction(total_potential, initial_conditions, 
                                                   0, timestep, N_timesteps);
    }

    std::vector<std::vector<vrs::Coords2d>>
    trajectories = multithreading::executeInParallel(tp_integration_functions);
    
    utility::writeCsv("../data/test_particle_1.csv", vrs::getPolarsFlat(utility::flatten(trajectories)));

}

    