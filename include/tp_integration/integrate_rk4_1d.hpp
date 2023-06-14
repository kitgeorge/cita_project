





namespace tp_integration {

std::function<std::vector<std::array<double, 2>>()>
getTpIntegrationFunction1d(const std::function<double(double)>& force,
                           std::array<double, 2> initial_conditions,
                           double timestep, double integration_time);

std::array<double, 2>
rk4Iteration1d(const std::function<double(double)>& force,
               std::array<double, 2> coords, double timestep);


}