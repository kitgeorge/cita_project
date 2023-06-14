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
#include "mapXVtoAA2D.hpp"
#include "dehnen_df.hpp"
#include "theta_r_integrator.hpp"
#include <random>

namespace vrs = vectors;
namespace ptl = potential;

int main() {
    int N_test_particles = 1e3;
    int N_timesteps = 1e4;

    double mestel_vc = 220*Units::kms;
    double R_0 = 8*Units::kpc;
    double spiral_amplitude_fraction = 0;
    double local_pitch_angle = std::numbers::pi/6;
    int m = 2;
    double k_R = m/R_0/tan(local_pitch_angle);
    double local_Omega = mestel_vc/R_0;
    double local_kappa = local_Omega*sqrt(2);
    double pattern_speed = local_Omega + local_kappa/m;

    double integration_time = 50*2*std::numbers::pi/local_Omega;
    double timestep = integration_time/N_timesteps;

    ptl::AxsymFuncs axsym_potential = ptl::getMestel(mestel_vc, R_0);
    ptl::PotentialFuncs total_potential = axsym_potential;
    int N_spirals = 1;

    // std::cout << total_potential.potential(R_0, 0, 0.05*integration_time) << ", "
    //           << total_potential.polar_force(R_0, 0, 0.05*integration_time)[0] << ", "
    //           << total_potential.cartesian_force(R_0, 0, 0.05*integration_time)[0]
    //           << std::endl; 

    for(int i = 0; i < N_spirals; ++i) {
        // int i = 0;
        double phase = (double)std::rand()/RAND_MAX*2*std::numbers::pi;
        ptl::PotentialFuncs spiral_potential 
        = ptl::getSpiralPotential(m, k_R, spiral_amplitude_fraction/k_R
                                            *(-1)*axsym_potential.polar_force(R_0, 0, 0)[0],
                                pattern_speed, phase);
        // std::function<double(double, double, double)>
        // envelope = [i, N_spirals, integration_time] (double R, double phi, double t) {
        //     double t_max = ((double)i + 0.5)/N_spirals*integration_time;
        //     double sigma_t = integration_time/(2*N_spirals);
        //     return exp(-pow(t - t_max, 2)/(2*pow(sigma_t, 2)));
        // };

        // spiral_potential.multiply(envelope);

        // ptl::PotentialFuncs total_potential(axsym_potential, spiral_potential);
        total_potential += spiral_potential;
    }

    ////////////
    // Debugging code
    //////////// 

    // int N_R_samples = 100;
    // int N_phi_samples = 100;
    // int N_t_samples = 100;
    // double R_max = 2*R_0;

    // std::vector<std::vector<std::vector<std::array<double, 2>>>>
    // force_difference(N_R_samples,
    //     std::vector<std::vector<std::array<double, 2>>>(N_phi_samples,
    //         std::vector<std::array<double, 2>>(N_t_samples)));
    // for(int i = 0; i < N_R_samples; ++i) {
    //     for(int j = 0; j < N_phi_samples; ++j) {
    //         for(int k = 0; k < N_t_samples; ++k) {
    //             double R = (double)(i + 1)/N_R_samples*R_max;
    //             double phi = (double)j/N_phi_samples*2*std::numbers::pi;
    //             double t = (double)k/N_t_samples*integration_time;
    //             // force_difference[i][j][k] = {{total_potential.polar_force(R, phi, t)[0] - axsym_potential.polar_force(R, phi, t)[0],
    //                                         //   total_potential.polar_force(R, phi, t)[1] - axsym_potential.polar_force(R, phi, t)[1]}};
    //             force_difference[i][j][k] = {{total_potential.cartesian_force(R*cos(phi), R*sin(phi), t)[0] - axsym_potential.cartesian_force(R*cos(phi), R*sin(phi), t)[0],
    //                                           total_potential.cartesian_force(R*cos(phi), R*sin(phi), t)[1] - axsym_potential.cartesian_force(R*cos(phi), R*sin(phi), t)[1]}};
    //         }
    //     }
    // }

    // utility::writeCsv("../data/debug.csv", 
    //                   utility::flatten(utility::flatten(force_difference)));
    // return 0;
    //////////////
    //////////////


    double gamma = 2*local_Omega/local_kappa;

    // Profiles to give a constant Q
    std::function<double(double)>
    surface_density = [R_0] (double R) {
        return 49*Units::Msun/pow(Units::pc, 2)*exp(-R/R_0);
    };
    std::function<double(double)>
    sigma_R = [surface_density, axsym_potential] (double R) {
        double Q = 2;
        return 3.36*Units::G*Q*surface_density(R)/axsym_potential.kappa(R);
    };

    // std::array<std::array<double, 2>, 2> 
    // bounds = {{ {{axsym_potential.potential_R(0.5*R_0),
    //               axsym_potential.potential_R(2*R_0)}},
    //             {{axsym_potential.LcGivenRc(0.5*R_0),
    //               axsym_potential.LcGivenRc(2*R_0)}} }};

    std::array<std::array<double, 2>, 2> 
    bounds = {{ {{axsym_potential.EcGivenRc(0.9*R_0),
                  axsym_potential.EcGivenRc(1.1*R_0)}},
                       {{axsym_potential.LcGivenRc(0.9*R_0),
                         axsym_potential.LcGivenRc(1.1*R_0)}} }};             

    std::function<double(double, double)>
    dehnen_df = df::getNormDehnenDF(axsym_potential, surface_density,
                                    sigma_R, bounds[0]);
    
    std::vector<vrs::Coords2d>
    initial_conditions(N_test_particles);

    double u_max = 5;
    int N_intervals = 1e3;
    int N_iterate = 5;

    for(int i = 0; i < N_test_particles; ++i) {
        std::array<double, 2> 
        sample_E_L = df::getDFSampleEL(dehnen_df, bounds);
        double theta_R = (double)std::rand()/RAND_MAX*2*std::numbers::pi;
        actions::ThetaRIntegrator integrator(axsym_potential, sample_E_L[0],
                                             sample_E_L[1], 
                                             u_max, N_intervals, N_iterate);
        std::array<std::array<double, 2>, 2> init_coords;
        std::array<double, 2> R_coords = integrator.getCoords(theta_R);
        init_coords[0][0] = R_coords[0];
        init_coords[1][0] = R_coords[1];
        init_coords[0][1] = (double)std::rand()/RAND_MAX*2*std::numbers::pi;
        init_coords[1][1] = sample_E_L[1]/init_coords[0][0];
        initial_conditions[i].setPolar(init_coords);
    }

    utility::writeCsv("../data/tp_dehnen_df_init.csv", vrs::getPolarsFlat(initial_conditions));



    std::vector<std::function<std::vector<vrs::Coords2d>()>>
    tp_integration_functions(N_test_particles);
    for(int i = 0; i < N_test_particles; ++i) {
        tp_integration_functions[i]
        = tp_integration::getTpIntegrationFunction(total_potential, initial_conditions[i],
                                                   0, timestep, N_timesteps);
    }

    std::vector<std::vector<vrs::Coords2d>>
    trajectories = multithreading::executeInParallel(tp_integration_functions);


    int N_timesteps_out = 500;
    std::vector<std::vector<vrs::Coords2d>>
    output_trajectories(N_test_particles,
                        std::vector<vrs::Coords2d>(N_timesteps_out));
    for(int i = 0; i < N_test_particles; ++i) {
        for(int j = 0; j < N_timesteps_out; ++j) {
            int k = j*N_timesteps/N_timesteps_out;
            output_trajectories[i][j] = trajectories[i][k];
        }
    }

    utility::writeCsv("../data/test_particle_dehnen.csv", vrs::getPolarsFlat(utility::flatten(output_trajectories)));
    return 0;
    RC::MapXVtoAA2D action_map(&axsym_potential);
    int N_tau = 1e3;
    std::vector<std::vector<std::vector<std::vector<double>>>>
    aa_trajectories(N_test_particles, 
                    std::vector<std::vector<std::vector<double>>>(N_timesteps_out,
                    std::vector<std::vector<double>>(2, std::vector<double>(2))));
    for(int i = 0; i < N_test_particles; ++i) {
        for(int k = 0; k < N_timesteps_out; ++k) {
            int j = k*N_timesteps/N_timesteps_out;
            std::cout << i << ", " << k << std::endl;
            action_map.mapXVtoAA(trajectories[i][j].polar[0][0],
                                 trajectories[i][j].polar[0][1],
                                 trajectories[i][j].polar[1][0],
                                 trajectories[i][j].polar[1][1],
                                 N_tau);
            aa_trajectories[i][k][0][0] = action_map.Jr;
            aa_trajectories[i][k][0][1] = action_map.Jpsi;
            aa_trajectories[i][k][1][0] = action_map.thetar;
            aa_trajectories[i][k][1][1] = action_map.thetapsi;
        }
    }

    utility::writeCsv("../data/test_particle_dehnen_aa.csv", utility::flatten(aa_trajectories));


}