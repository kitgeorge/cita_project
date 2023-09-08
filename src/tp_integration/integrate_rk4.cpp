#include "integrate_rk4.hpp"
#include <iostream>
#include <mutex>
#include <ctime>

namespace vrs = vectors;
namespace ptl = potential;

std::mutex rk4_mtx;

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
        assert(std::isfinite(output[0].polar[0][0]));
        double t = t_start;
        for(int i = 0; i < N_timesteps; ++i) {
            output[i + 1] = rk4Iteration(potential, output[i], t, timestep);
            t += timestep;
        }
        return output;
    };
}

std::vector<std::function<std::vector<vrs::Coords2d>()>>
getTpIntegrationFunctions(const potential::PotentialFuncs& potential,
                          const std::vector<vectors::Coords2d>& initial_conditions,
                          double t_start, double timestep, int N_timesteps) {
    int N_test_particles = initial_conditions.size();
    std::vector<std::function<std::vector<vrs::Coords2d>()>>
    output(N_test_particles);
    for(int i = 0; i < N_test_particles; ++i) {
        output[i] = getTpIntegrationFunction(potential, initial_conditions[i],
                                             t_start, timestep, N_timesteps);
    }
    return output;
}

vrs::Coords2d
rk4IterationBoxed(const ptl::PotentialFuncs& potential,
                  vrs::Coords2d coords, double t, double timestep,
                  double R_Ka) {
    double time_0 = std::time(nullptr);
    assert(std::isfinite(coords.polar[0][0]));
    assert(coords.polar[0][0] < 20*Units::kpc); // For debugging
    
    std::array<std::array<double, 2>, 2>
    cart = coords.cartesian;

    rk4_mtx.lock();
    // std::cout << "rk4Iteration, " << coords.polar[0][0] << ", " 
    //           << coords.polar[0][1] << ", " << cart[0][0] << ", "
    //           << cart[0][1] << std::endl;
    // std::cout << "velocities, " << coords.polar[1][0] << ", " 
    //           << coords.polar[1][1] << ", " << cart[1][0] << ", " 
    //           << cart[1][1] << std::endl;
    assert(std::isfinite(cart[0][0]));
    assert(std::isfinite(cart[0][1]));
    assert(std::isfinite(cart[1][0]));
    assert(std::isfinite(cart[1][1]));
    rk4_mtx.unlock();

    double x_temp;
    double y_temp;

    std::function<vrs::Coords2d()>
    reflect = [&coords]() {
        std::array<std::array<double, 2>, 2>
        polar = {{coords.polar[0], {{-coords.polar[1][0],
                                     coords.polar[1][1]}}}};
        return vrs::Coords2d(polar, 1);
    };

    std::array<std::array<std::array<double, 2>, 2>, 4> k;
    k[0][0] = cart[1];
    double time_1 = std::time(nullptr);
    k[0][1] = potential.cartesian_force(cart[0][0], cart[0][1], t);
    double time_2 = std::time(nullptr);
    rk4_mtx.lock();
    std::cout << "First force: " << time_2 - time_1 << std::endl;
    rk4_mtx.unlock();
    k[1][0] = add_arrays(cart[1], multiply_array(k[0][1], timestep/2));
    x_temp = cart[0][0] + timestep/2*k[0][0][0];
    y_temp = cart[0][1] + timestep/2*k[0][0][1];
    if(pow(x_temp, 2) + pow(y_temp, 2) > pow(R_Ka, 2)) {
        return reflect();
    }
    time_1 = std::time(nullptr);
    k[1][1] = potential.cartesian_force(x_temp, y_temp, t + timestep/2);
    time_2 = std::time(nullptr);
    rk4_mtx.lock();
    std::cout << "Second force: " << time_2 - time_1 << std::endl;
    rk4_mtx.unlock();
    k[2][0] = add_arrays(cart[1], multiply_array(k[1][1], timestep/2));
    x_temp = cart[0][0] + timestep/2*k[1][0][0];
    y_temp = cart[0][1] + timestep/2*k[1][0][1];
    if(pow(x_temp, 2) + pow(y_temp, 2) > pow(R_Ka, 2)) {
        return reflect();
    }
    time_1 = std::time(nullptr);
    k[2][1] = potential.cartesian_force(x_temp,
                                        y_temp,
                                        t + timestep/2);
    time_2 = std::time(nullptr);
    rk4_mtx.lock();
    std::cout << "Third force: " << time_2 - time_1 << std::endl;
    rk4_mtx.unlock();
    k[3][0] = add_arrays(cart[1], multiply_array(k[2][1], timestep));
    x_temp = cart[0][0] + timestep*k[2][0][0];
    y_temp = cart[0][1] + timestep*k[2][0][1];
    if(pow(x_temp, 2) + pow(y_temp, 2) > pow(R_Ka, 2)) {
        return reflect();
    }
    time_1 = std::time(nullptr);
    k[3][1] = potential.cartesian_force(x_temp, y_temp,
                                        t + timestep);
    time_2 = std::time(nullptr);
    rk4_mtx.lock();
    std::cout << "Fourth force: " << time_2 - time_1 << std::endl;
    rk4_mtx.unlock();

    cart = add_arrays(cart, multiply_array(k[0], timestep/6));
    cart = add_arrays(cart, multiply_array(k[1], timestep/3));
    cart = add_arrays(cart, multiply_array(k[2], timestep/3));
    cart = add_arrays(cart, multiply_array(k[3], timestep/6));

    vrs::Coords2d output(cart, 0);
    rk4_mtx.lock();
    std::cout << "Total time: " << std::time(nullptr) - time_0 << std::endl;
    rk4_mtx.unlock();
    return output;
}

vrs::Coords2d
rk4Iteration(const ptl::PotentialFuncs& potential, 
              vrs::Coords2d coords,
              double t, double timestep) {

    assert(std::isfinite(coords.polar[0][0]));
    assert(coords.polar[0][0] < 20*Units::kpc); // For debugging
    
    std::array<std::array<double, 2>, 2>
    cart = coords.cartesian;

    rk4_mtx.lock();
    // std::cout << "rk4Iteration, " << coords.polar[0][0] << ", " 
    //           << coords.polar[0][1] << ", " << cart[0][0] << ", "
    //           << cart[0][1] << std::endl;
    // std::cout << "velocities, " << coords.polar[1][0] << ", " 
    //           << coords.polar[1][1] << ", " << cart[1][0] << ", " 
    //           << cart[1][1] << std::endl;
    assert(std::isfinite(cart[0][0]));
    assert(std::isfinite(cart[0][1]));
    assert(std::isfinite(cart[1][0]));
    assert(std::isfinite(cart[1][1]));
    rk4_mtx.unlock();

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