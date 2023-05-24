#include <numbers>
#include "flatten.hpp"
#include "vector_io.hpp"
#include "add_functions.hpp"
#include "units.hpp"


int main() {
    int N_test_particles = 1e6;
    int N_timesteps = 1e3;
    double integration_time = 10;
    double timestep = integration_time/N_timesteps;

    double mestel_vc = 220*consts::km;
    double R_0 = 8*consts::kpc;
    double spiral_amplitude_fraction = 0.1;
    double local_pitch_angle = std::numbers::pi/6;
    int m = 2;
    double k_R = m/R_0/tan(local_pitch_angle);
    double local_Omega = mestel_vc/R_0;
    double local_kappa = local_Omega*sqrt(2);
    double pattern_speed = local_Omega + local_kappa/m;

    PotentialFuncs axsym_potential = getMestel(mestel_vc);
    PotentialFuncs spiral_potential 
    = getSpiralPotential(m, k_R, spiral_amplitude_fraction
                                 *axsym_potential.potential(R_0),
                         pattern_speed);
    PotentialFuncs total_potential(axsym_potential, spiral_potential);
    

    double sigma_R = 30*consts::km;
    double gamma = 2*local_Omega/local_kappa;
    std::function<double(double, double)>
    schwarzchild_df = [=] (double v_R, double v_phi) {
        return exp(-(pow(v_R, 2) + pow(gamma*(v_phi - v_c), 2))/(2*pow(sigma_R, 2)));
    };
    double R_min = 0.9*R_0;
    double R_max = 1.1*R_0;
    double v_R_max = 100*consts::km;
    double v_phi_max = v_R_max/gamma;








    std::vector<std::function<std::vector<std::array<std::array<double, 2>, 2>>()>
    tp_integration_functions(N_test_particles); 

    for(int i = 0; i < N_test_particles; ++i
        std::array<std::array<double, 2>, 2> 
        initial_conditions = getDFSample(schwarzchild_df, 
                                         R_min, R_max, 
                                         v_R_max, v_phi_max);
        tp_integration_functions[i] 
        = getTpIntegrationFunction(potential, initial_conditions, N_timesteps);
    }

    std::vector<std::vector<std::array<std::array<double, 2>, 2>>>
    trajectories = multithreading::execute_in_parallel(tp_integration_functions);
    
    writeCsv("test_particle_1.csv", utility::flatten(trajectories));

}

    