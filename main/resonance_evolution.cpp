#include "dehnen_df.hpp"
#include "libration_calculator.hpp"
#include "mestel.hpp"
#include "mestel_spiral.hpp"
#include "units.hpp"
#include "mapXVtoAA2D.hpp"
#include "vector_io.hpp"

int main() {


    // potential::MestelObjects mest = potential::getStandardSpiral();
    potential::MestelObjects mest = potential::getGlobalSpiral();
    // OLR
    int l = 1;
    int m = mest.N_phi;
    double R_0 = mest.R_res[2];
    actions::JresFinder finder = mest.res_finders[2].value();

    std::function<std::array<double, 2>(std::array<double, 2>)>
    J_s_f = [l, m] (std::array<double, 2> actions) {
        // {J_s, J_f as a function of J_phi, J_R}
        std::array<double, 2> output = {actions[0]/m, actions[1] - actions[0]*l/m}; 
        return output;
    };

    double J_phi_0 = m*finder.J_s_res_0;
    std::array<std::array<double, 2>, 2>
    bounds = {{ {{0.9*J_phi_0, 1.05*J_phi_0}},
                {{0, 0.05*J_phi_0}} }};
    double dJ = J_phi_0/1000;
    int N_J_phi = (bounds[0][1] - bounds[0][0])/dJ;
    int N_J_R = (bounds[1][1] - bounds[1][0])/dJ;

    double delta_J_s = 0.001*finder.J_s_res_0;
    int N_jacobi_intervals = 1e5;
    
////////////////////////////////////////////
////// Dehnen DF (in (E, J) and in (J_phi, J_R))
////////////////////////////////////////////

    // Profiles to give a constant Q
    std::function<double(double)>
    surface_density = [R_0] (double R) {
        return 49*Units::Msun/pow(Units::pc, 2)*exp(-R/R_0);
    };
    potential::AxsymFuncs axsym_potential = mest.mestel;
    std::function<double(double)>
    sigma_R = [surface_density, axsym_potential] (double R) {
        double Q = 2;
        return 3.36*Units::G*Q*surface_density(R)/axsym_potential.kappa(R);
    };
    std::function<double(double, double)>
    initial_df = df::getDehnenDF(axsym_potential, surface_density, sigma_R);

    bool prog = 1;
    int N_tau = 1000;
    std::vector<std::vector<double>> E_J;
    RC::gridEOverJ(&mest.mestel, N_tau, (bounds[1][1])/dJ, 
                   (bounds[0][1])/dJ, dJ, dJ, E_J, prog);
    std::function<double(double, double)>
    initial_df_actions = [E_J, initial_df, dJ] (double J_phi, double J_R) {
        double E = RC::XGivenJ(J_R, J_phi, dJ, dJ, E_J);
        return initial_df(E, J_phi);
    };

////////////////////////////////////////////
////// Calculate time-evolution of librating point in phase space
////////////////////////////////////////////
    
    int N_angle_samples = 10;

    std::function<double(std::array<double, 2>, double)>
    get_evolved_df = [&](std::array<double, 2> actions,
                         double t) {
        // For now, we'll only consider libration and not circulation
        std::array<double, 2>
        slow_fast = J_s_f(actions);
        actions::LibrationCalculator 
        calc(slow_fast[1], finder, mest.libfuncs, delta_J_s, N_jacobi_intervals);
        std::vector<double> slow_angles(N_angle_samples);
        for(int i = 0; i < N_angle_samples; ++i) {
            slow_angles[i] = -calc.g -std::numbers::pi
                             + (double)i/N_angle_samples*2*std::numbers::pi;
        }
        double output = 0;
        for(int i = 0; i < N_angle_samples; ++i) {
            double J_s_old;
            if(slow_fast[0] > calc.separatrices(slow_angles[i])[0]
               && slow_fast[0] < calc.separatrices(slow_angles[i])[1]) {
                calc.setICs(slow_fast[0], slow_angles[i]);
                J_s_old = calc.trajectory(-t)[1];
            }
            else {
                J_s_old = slow_fast[0];
            }
            double J_phi_old = m*J_s_old;
            double J_R_old = l*J_s_old + slow_fast[1];
            output += initial_df_actions(J_phi_old, J_R_old)/N_angle_samples;
        }
        return output;
    };

//////////////////////////////////////////////
////// Constructing functions for each (J_phi, J_R, t) (to execute in parallel)
//////////////////////////////////////////////

    double timestep = 50*Units::Myr;
    int N_timesteps = 10;

    std::vector<std::function<double()>>
    evolved_df_functions(N_J_phi*N_J_R*N_timesteps);
    for(int i = 0; i < N_timesteps; ++i) {
        for(int j = 0; j < N_J_phi; ++j) {
            for(int k = 0; k < N_J_R; ++k) {
                double t = i*timestep;
                double J_phi = bounds[0][0] + j*dJ;
                double J_R = bounds[1][0] = k*dJ;
                evolved_df_functions[i*N_J_phi*N_J_R + j*N_J_R + k]
                = [t, J_phi, J_R, get_evolved_df, i, j, k] () {
                    std::cout << i << ", " << j << ", " << k << std::endl;
                    return get_evolved_df({{J_phi, J_R}}, t);
                };
            }
        }
    }

/////////////////////////////////////////////
/////////////////////////////////////////////

    std::vector<double> 
    evolved_df_values = multithreading::executeInParallel(evolved_df_functions);

    utility::writeCsv("../data/resonance_evolution_global.csv", evolved_df_values);

}