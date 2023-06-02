#include "integrate_rk4.hpp"

namespace vrs = vectors;
namespace ptl = potential;

namespace tp_integration {

using utility::add_arrays;
using utility::multiply_array;

std::function<std::vector<vrs::Coords2d>()>
getTpIntegrationFunction(const ptl::PotentialFuncs& potential,
                         const vrs::Coords2d& initial_conditions,
                         double t_start, 
                         double timestep, int N_timesteps) {
    return [=] () {
        std::vector<vrs::Coords2d>
        output(N_timesteps + 1);
        output[0] = initial_conditions;
        double t = t_start;
        for(int i = 0; i < N_timesteps; ++i) {
            output[i + 1] = rk4Iteration(potential, output[i], t, timestep);
            t += timestep;
        }
        return output;
    };
}

vrs::Coords2d
rk4Iteration(const ptl::PotentialFuncs& potential, 
              vrs::Coords2d coords,
              double t, double timestep) {
    std::array<std::array<double, 2>, 2>
    cart = coords.cartesian;

    std::array<std::array<std::array<double, 2>, 2>, 4> k;
    k[0][0] = cart[1];
    k[0][1] = potential.cartesian_force(cart[0][0], cart[0][1], t);
    k[1][0] = add_arrays(cart[1], multiply_array(k[0][1], timestep/2));
    k[1][1] = potential.cartesian_force(cart[0][0] + timestep/2*k[0][0][0],
                                        cart[0][1] + timestep/2*k[0][0][1],
                                        t + timestep/2);
    k[2][0] = add_arrays(cart[1], multiply_array(k[1][1], timestep/2));
    k[2][1] = potential.cartesian_force(cart[0][0] + timestep/2*k[1][0][0],
                                        cart[0][1] + timestep/2*k[1][0][1],
                                        t + timestep/2);
    k[3][0] = add_arrays(cart[1], multiply_array(k[2][1], timestep));
    k[3][1] = potential.cartesian_force(cart[0][0] + timestep*k[2][0][0],
                                        cart[0][1] + timestep*k[2][0][1],
                                        t + timestep);

    cart = add_arrays(cart, multiply_array(k[0], timestep/6));
    cart = add_arrays(cart, multiply_array(k[1], timestep/3));
    cart = add_arrays(cart, multiply_array(k[2], timestep/3));
    cart = add_arrays(cart, multiply_array(k[3], timestep/6));

    vrs::Coords2d output(cart, 0);

    return output;
}

}