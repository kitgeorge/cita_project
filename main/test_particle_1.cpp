int main() {
    std::function<double(double, double)>
    axsym_potential = getAxsymPotential(axsym_parameters);
    std::function<double(double, double)>
    spiral_potential = getSpiralPotential(spiral_parameters);
    std::function<double(double, double)>
    potential = add_functions(axsym_parameters, spiral_parameters);
    
    std::vector<std::function<std::vector<std::array<double, 2>>()>
    tp_integration_functions(N_test_particles); 

    for(int i = 0; i < N_test_particles; ++i) {
        std::array<double, 2> 
        initial_conditions = getInitialConditions(ic_params, i);
        tp_integration_functions[i] 
        = getTpIntegrationFunction(potential, initial_conditions, N_timesteps);
    }

    std::vector<std::vector<std::array<double, 2>>>
    trajectories = execute_in_parallel(tp_integration_functions);
    
    writeCsv("test_particle_1.csv", flatten(trajectories));

}

    