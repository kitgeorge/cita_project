#include "add_arrays.hpp"
#include "coords.hpp"
#include "potential_funcs.hpp"



namespace tp_integration {

std::function<std::vector<vrs::Coords2d>()>
getTpIntegrationFunction(const potential::PotentialFuncs& potential,
                         const vectors::Coords2d& initial_conditions,
                         double t_start, 
                         double timestep, int N_timesteps);

std::vector<std::function<std::vector<vrs::Coords2d>()>>
getTpIntegrationFunctions(const potential::PotentialFuncs& potential,
                          const std::vector<vectors::Coords2d>& initial_conditions,
                          double t_start, double timestep, int N_timesteps);

vectors::Coords2d
rk4Iteration(const potential::PotentialFuncs& potential, 
              vectors::Coords2d coords,
              double t, double timestep);

vrs::Coords2d
rk4IterationBoxed(const potential::PotentialFuncs& potential,
                  vrs::Coords2d coords, double t, double timestep,
                  double R_Ka);





}