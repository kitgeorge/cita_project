#include "bfe_nbody.hpp"

namespace n_body {

std::mutex mtx;

void BFENBody::getPotential() {
    std::vector<std::array<double, 3>> density(N_particles);
    for(int i = 0; i < N_particles; ++i) {
        density[i] = {coords[i][0][0], coords[i][0][1], masses[i]};
    }
    bfe_pot.initFromDensity(density);
}

potential::PotentialFuncs BFENBody::getInit() {
    getPotential();
    return potential::PotentialFuncs(bfe_pot);
}

void BFENBody::iterate() {
    potential::PotentialFuncs pot({background, init*(-1), 
                                    potential::PotentialFuncs(bfe_pot)});
    // Ensure that N_particles is a multiple of particles_per_function please
    int particles_per_function = 2;
    int N_functions = N_particles/particles_per_function;
    std::vector<std::function<std::vector<std::array<std::array<double, 2>, 2>>()>>
    rk4_iteration_functions(N_functions);
    for(int i = 0; i < N_functions; ++i) {
        std::cout << "Making RK4 function " << i << std::endl;
        rk4_iteration_functions[i] = [i, particles_per_function, &pot, coords=coords, timestep=timestep] () {
            mtx.lock();
            std::cout << "Executing RK4 function " << i << std::endl;
            mtx.unlock();
            std::vector<std::array<std::array<double, 2>, 2>>
            output(particles_per_function);
            for(int j = 0; j < particles_per_function; ++j) {
                // Should link this to BFE R_Ka somehow
                double R_Ka = 20*Units::kpc;
                std::array<std::array<double, 2>, 2>
                coords_ = coords[i*particles_per_function + j];
                // Keep particles within R_Ka, approximately by
                // reflecting them off a wall at R_Ka
                if(coords_[0][0] > R_Ka) {
                    coords_[0][0] = 0.999*R_Ka;
                    coords_[1][0] = -std::abs(coords_[1][0]);
                }
                vectors::Coords2d 
                x(coords_, 1);
                assert(std::isfinite(x.polar[0][0]));
                assert(x.polar[0][0] < R_Ka);
                output[j] = tp_integration::
                            rk4IterationBoxed(pot, x, 0, timestep, R_Ka).polar;
            }
            return output;
        };
    }
    std::vector<std::vector<std::array<std::array<double, 2>, 2>>>
    data = multithreading::executeInParallel(rk4_iteration_functions);
    coords = utility::flatten(data);
    getPotential();
}

BFENBody::BFENBody(double timestep_, int save_interval_, 
            double integration_time, int N_particles_,
            const potential::AxsymFuncs& background_,
            const std::vector<double>& masses_,
            const std::vector<std::array<std::array<double, 2>, 2>>&
            init_coords): expansion(std::make_shared<const basis_functions::BFE>()),
            timestep(timestep_), 
            save_interval(save_interval_),
            N_timesteps(integration_time/timestep_), 
            N_particles(N_particles_), background(background_),
            masses(masses_), 
            coords(init_coords),
            bfe_pot(expansion),
            init(getInit()),
            bfe_coefficients(N_timesteps + 1),
            bfe_coefficient_norms(N_timesteps + 1) {
    std::array<int, 2> shape = {{N_timesteps/save_interval + 1,
                                                N_particles}};
    std::cout << "Setting up trajectories table..." << std::endl;
    saved_trajectories = utility::
        makeShape<std::array<std::array<double, 2>, 2>>(shape);
    std::cout << "Tabulating initial coords in trajectories table..." 
              << std::endl;
    saved_trajectories[0] = coords;
    std::cout << "Tabulating initial bfe coefficients" << std::endl;
    bfe_coefficients[0] = bfe_pot.getCoefficients();
    std::cout << "Tabulating initial bfe coefficient norms" << std::endl;
    bfe_coefficient_norms[0] = bfe_pot.calculateAbsNorm();
    for(int i = 0; i < N_timesteps; ++i) {
        std::cout << "Iteration " << i << std::endl;
        iterate();
        if(i % save_interval == 0) {
            saved_trajectories[i/save_interval] = coords;
        }
        bfe_coefficients[i] = bfe_pot.getCoefficients();
        bfe_coefficient_norms[i] = bfe_pot.calculateAbsNorm();
    }
}

utility::vector2d<std::array<std::array<double, 2>, 2>>
BFENBody::getTrajectories() {
    return saved_trajectories;
}

utility::vector3d<std::complex<double>>
BFENBody::getBFECoefficients() {
    return bfe_coefficients;
}

std::vector<double>
BFENBody::getBFECoefficientNorms() {
    return bfe_coefficient_norms;
}

}