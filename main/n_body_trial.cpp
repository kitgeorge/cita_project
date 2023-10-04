#define NDEBUG
#include "bfe_nbody.hpp"
#include "units.hpp"
#include "tapered_df.hpp"
#include "mestel.hpp"
#include "sample.hpp"
#include "vector_io.hpp"


int main() {
    utility::open_channel(1);
    double timestep = 1*Units::Myr;
    int save_interval = 1;
    double integration_time = 1000*Units::Myr;

    double v_c = 220*Units::kms;
    double R_0 = 1*Units::kpc;
    double active_fraction = 0.5;
    std::array<double, 2> taper_radii = {1, 10};
    std::array<double, 2> taper_indices = {4, 5};
    std::array<double, 2> cutoff_radii = {0.1, 20};
    std::function<double(double)>
    target_Q = [] (double R) {return 2;};
    df::TaperedDF tap(v_c, R_0, active_fraction,
                  taper_radii, taper_indices, cutoff_radii,
                  target_Q);
    std::array<std::array<double, 2>, 2>
    E_L_bounds = {tap.E_bounds, tap.L_bounds};

    double u_max = 3;
    int N_u_intervals = 1000;
    int N_u_iterate = 5;

    int N_particles = 1e4;
    
    std::vector<std::array<std::array<double, 2>, 2>>
    sample_coords(N_particles);
    std::vector<std::array<double, 2>> sample_positions(N_particles);

    for(int i = 0; i < N_particles; ++i) {
        // For some reason, a sample was giving R = -nan. Rather than
        // fixing the problem I'm just going to suppress it for now
        int nan_flag = 1;
        while(nan_flag) {
            sample_coords[i] = df::getDFSampleViaEL(tap.getTaperedDF(), E_L_bounds,
                                    potential::getMestel(v_c, R_0),
                                    u_max, N_u_intervals, N_u_iterate);
            if(std::isfinite(sample_coords[i][1][0])) {
                nan_flag = 0;
            }
        }
        sample_positions[i] = sample_coords[i][0];
    }

    ///////////////////////////////////////////////////////////////////
    /// Debugging sample angles
    ///////////////////////////////////////////////////////////////////

    // potential::AxsymFuncs pot = potential::getMestel(v_c, R_0);
    // std::vector<double> theta_R_values(N_particles);
    // for(int i = 0; i < N_particles; ++i) {
    //     double E = pot.EGivenPolar(sample_coords[i]);
    //     double L = pot.LGivenPolar(sample_coords[i]);
    //     actions::ThetaRIntegrator integrator(pot, E, L, u_max, N_u_intervals, N_u_iterate);
    //     theta_R_values[i] = integrator.calculateThetaR({{sample_coords[i][0][0],
    //                                                      sample_coords[i][1][0]}});
    // }

    // utility::writeCsv("../data/n_body/init_theta_R_values.csv",
    //                   theta_R_values);
    
    ///////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////

    

    potential::AxsymFuncs background = potential::getMestel(v_c, R_0);


    std::function<double(double, double)>
    target_density = [v_c, active_fraction] (double R, double phi) {
        return active_fraction*pow(v_c, 2)/(2*std::numbers::pi*Units::G*R);
    };

    double particle_mass = basis_functions::
        getParticleMass(sample_positions, target_density, {3*R_0, 0});

    std::vector<double> particle_masses(N_particles);
    for(int i = 0; i < N_particles; ++i) {
        particle_masses[i] = particle_mass;
    }

    n_body::BFENBody simulation(timestep, save_interval, integration_time,
                        N_particles, background, particle_masses, 
                        sample_coords);
    utility::writeCsv("../data/n_body/test_trajectories.csv",
                      utility::flatten(utility::flatten(utility::flatten(
                        simulation.getTrajectories()
                      ))));
    std::array<int, 3> 
    coefficients_shape = utility::getShape(simulation.getBFECoefficients());
    utility::vector3d<std::array<double, 2>>
    coefficients = utility::makeShape<std::array<double, 2>>(coefficients_shape);
    utility::vector3d<std::complex<double>> c_coefficients = simulation.getBFECoefficients();
    for(int i = 0; i < coefficients_shape[0]; ++i) {
        for(int j = 0; j < coefficients_shape[1]; ++j) {
            for(int k = 0; k < coefficients_shape[2]; ++k) {
                coefficients[i][j][k][0] = c_coefficients[i][j][k].real();
                coefficients[i][j][k][1] = c_coefficients[i][j][k].imag();
            }
        }
    }
    utility::writeCsv("../data/n_body/test_bfe_coefficients.csv",
                      utility::flatten(utility::flatten(coefficients)));
    utility::writeCsv("../data/n_body/test_bfe_coefficient_norms.csv",
                      simulation.getBFECoefficientNorms());

    ///////////////////////////////////////////////////////////////////
    /// Calculate angles throughout trajectories
    /// also E, L values
    ///////////////////////////////////////////////////////////////////

    potential::AxsymFuncs pot = potential::getMestel(v_c, R_0);
    int N_timesteps = integration_time/timestep;
    std::array<int, 2> shape = {N_particles, N_timesteps/save_interval + 1};
    // std::vector<std::function<std::vector<double>()>>
    // theta_R_functions(shape[0]);
    std::vector<std::function<std::vector<std::array<double, 2>>()>>
    E_L_functions(shape[0]);
    std::vector<std::vector<std::array<std::array<double, 2>, 2>>>
    R_coords_vector = simulation.getTrajectories();
    // std::vector<std::vector<std::array<double, 2>>>
    // E_L_values = utility::makeShape<std::array<double, 2>>(shape);
    for(int i = 0; i < N_particles; ++i) {
        E_L_functions[i] = [i, N_timesteps, save_interval, &pot, shape, &R_coords_vector] () {
            std::cout << "Calculating E, L: particle " << i << std::endl;
            std::vector<std::array<double, 2>> output(shape[1]);
            for(int j = 0; j <= N_timesteps/save_interval; ++j) {
                std::array<double, 2>
                output[j] = {{pot.EGivenPolar(R_coords_vector[j][i]),
                              pot.LGivenPolar(R_coords_vector[j][i])}};
                std::cout << output[j][0] << ", " << output[j][1] << std::endl;
            }
            return output;
        };
        // theta_R_functions[i] = [i, N_timesteps, save_interval, &pot, u_max, N_u_intervals, N_u_iterate, &simulation, shape, &R_coords_vector]() {
        //     std::cout << "Calculating angles: particle " << i << std::endl; 
        //     double E = pot.EGivenPolar(R_coords_vector[0][i]);
        //     double L = pot.LGivenPolar(R_coords_vector[0][i]);
        //     actions::ThetaRIntegrator
        //     integrator(pot, E, L, u_max, N_u_intervals, N_u_iterate);
        //     std::vector<double> output(shape[1]);
        //     for(int j = 0; j <= N_timesteps/save_interval; ++j) {
        //         std::array<double, 2>
        //         R_coords = {R_coords_vector[j][i][0][0],
        //                     R_coords_vector[j][i][1][0]};
        //         output[j] = integrator.calculateThetaR(R_coords);
        //     }
        //     return output;
        // };
    }
    utility::vector2d<std::array<double, 2>>
    E_L_values = multithreading::executeInParallel(E_L_functions);
    utility::writeCsv("../data/n_body/E_L_values.csv",
                      utility::flatten(utility::flatten(E_L_values)));

    ///////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////

}