#include "potential_from_density.hpp"
#include <iostream>
#include <mutex>
#include <chrono>

std::mutex pfd_mtx;

namespace basis_functions {


PotentialFromDensity::
PotentialFromDensity(int n_max, int l_max): 
                     expansion(std::make_shared<const basis_functions::BFE>()),
                     nl_max({n_max, l_max}) {}

PotentialFromDensity::
PotentialFromDensity(const std::shared_ptr<const basis_functions::BFE> expansion_,
                     int n_max, int l_max): expansion(expansion_),
                     nl_max({n_max, l_max}) {}

PotentialFromDensity::
PotentialFromDensity(const PotentialFromDensity& old):
                     expansion(old.expansion),
                     nl_max(old.nl_max), coefficients(old.coefficients),
                     trunc_density(old.trunc_density),
                     trunc_potential(old.trunc_potential),
                     trunc_force(old.trunc_force) {}

template <typename DensityType>
void PotentialFromDensity::
initFromDensity(const DensityType& density) {
    coefficients = getCoefficients(density);
    trunc_density = getTruncDensity();
    trunc_potential = getTruncPotential();
    trunc_force = getTruncForce();
}

void PotentialFromDensity::
initFromCoefficients(const std::vector<std::vector<std::complex<double>>>
                     coefficients_) {
    coefficients = coefficients_;
    trunc_density = getTruncDensity();
    trunc_potential = getTruncPotential();
    trunc_force = getTruncForce();
}


template <typename DensityType>
std::vector<std::vector<std::complex<double>>> 
PotentialFromDensity::
getCoefficients(const DensityType& density) const {
    pfd_mtx.lock();
    std::cout << nl_max[0] + 1 << ", " << nl_max[1] + 1 << std::endl;
    pfd_mtx.unlock();
    std::vector<std::vector<std::complex<double>>>
    output(nl_max[0] + 1, std::vector<std::complex<double>>(nl_max[1] + 1));
    std::vector<std::function<std::complex<double>()>>
    coefficient_functions((nl_max[0] + 1)*(nl_max[1] + 1));
    for(int i = 0; i <= nl_max[0]; ++i) {
        for(int j = 0; j <= nl_max[1]; ++j) {
            coefficient_functions[i*(nl_max[1] + 1) + j] = [i, j, expansion=expansion, &density, &pfd_mtx] () {
                // pfd_mtx.lock();
                // std::cout << "Calculating BFE coefficients: "
                //           << i << ", " << j << std::endl;
                // pfd_mtx.unlock();
                return expansion->getCoefficient(i, j, density);
            };
        }
    }
    std::vector<std::complex<double>>
    output_flat = multithreading::executeInParallel(coefficient_functions);
    for(int i = 0; i <= nl_max[0]; ++i) {
        for(int j = 0; j <= nl_max[1]; ++j) {
            output[i][j] = output_flat[i*(nl_max[1] + 1) + j];
        }
    }
    return output;
}

template <typename DataType>
std::vector<std::vector<std::function<DataType(double, double)>>>
PotentialFromDensity::
getTerms(std::function<std::function<DataType(double, double)>
                         (int, int)> BFE_member_function) const {
    auto time_0 = std::chrono::steady_clock::now();
    // Initialising table for 0 <= n <= n_max, -l_max <= 0 <= l_max
    std::vector<std::vector<std::function<DataType(double, double)>>>
    output(nl_max[0] + 1, 
            std::vector<std::function<DataType(double, double)>>(2*nl_max[1] + 1));
    // Populates table with respective basis functions multiplied by coefficients,
    // and complex conjugate for negative l (not sure if it should be a pointer
    // to expansion)
    for(int i = 0; i <= nl_max[0]; ++i) {
        for(int j = 0; j <= nl_max[1]; ++j) {
            // auto time_1 = std::chrono::steady_clock::now();
            output[i][nl_max[1] + j] = utility::multiplyFunction(BFE_member_function(i, j),
                                            coefficients[i][j]);
            // auto time_2 = std::chrono::steady_clock::now();
            // pfd_mtx.lock();
            // std::cout << "BFE terms: " << i << ", " << j << ", "
            //           << std::chrono::duration_cast<std::chrono::microseconds>
            //                     (time_2 - time_1).count() << "us"
            //           << std::endl;
            // pfd_mtx.unlock();
            if(j != 0) {
                output[i][nl_max[1] - j] 
                        = utility::conjugateFunction(output[i][nl_max[1] + j]);
            }
        }
    }
    // pfd_mtx.lock();
    // std::cout << "Total time: " 
    //           << std::chrono::duration_cast<std::chrono::milliseconds>
    //                             (std::chrono::steady_clock::now() - time_0).count()
    //           << "ms" << std::endl;
    // pfd_mtx.unlock();
    return output;
}

// Trying to generalise the summation of basis functions for density, potential and
// force (in some ways this is quite ugly. DataType refers to complex or 
// complex array (force))
template <typename DataType>
std::function<DataType(double, double)>
PotentialFromDensity::
getTruncFunction(std::function<std::function<DataType(double, double)>
                    (int, int)> BFE_member_function) const {
    // auto output_function = utility::addFunctions(
    //             utility::flatten(getTerms(BFE_member_function)));

    // auto debug_output = [output_function] (double R, double phi) {
    //     double x = 43543.554543;
    //     auto time_a = std::chrono::steady_clock::now();
    //     double y = x;
    //     auto time_b = std::chrono::steady_clock::now();
    //     // std::cout << "Time to read double: "
    //     //         << std::chrono::duration_cast<std::chrono::nanoseconds>
    //     //                 (time_b - time_a).count() << "ns"
    //     //         << std::endl;
    //     auto time_0 = std::chrono::steady_clock::now();
    //     auto output = output_function(R, phi);
    //     auto time_1 = std::chrono::steady_clock::now();
    //     // pfd_mtx.lock();
    //     // std::cout << "Total force time: " 
    //     //         << std::chrono::duration_cast<std::chrono::milliseconds>
    //     //                 (time_1 - time_0).count() << "ms"
    //     //         << std::endl;
    //     // pfd_mtx.unlock();
    //     return output;
    // };
    // return debug_output;

    // Would just do nested for loops for n and l, but we need to initialise output
    // with the first function before we can pass it to addFunctions.
    std::function<DataType(double, double)> 
    output = utility::multiplyFunction(BFE_member_function(0, 0), coefficients[0][0]);
    for(int l = 1; l <= nl_max[1]; ++l) {
        output = utility::addFunctions(output,
                                        utility::multiplyFunction(BFE_member_function(0, l),
                                            coefficients[0][l]));
    }
    for(int n = 1; n <= nl_max[0]; ++n) {
        for(int l = 0; l <= nl_max[1]; ++l) {
            output = utility::addFunctions(output,
                                           utility::multiplyFunction(BFE_member_function(n, l),
                                                coefficients[n][l]));
        }
    }

    return output;


    // Usually this is all
    // return utility::addFunctions(
                // utility::flatten(getTerms(BFE_member_function)));

}

std::function<double(double, double)>
PotentialFromDensity::
getTruncDensity() const {
    std::function<std::function<std::complex<double>(double, double)>(int, int)>
    wrapper = [=] (int i, int j) {
        return expansion->rho(i, j);
    };
    return utility::realFunction(getTruncFunction(wrapper));
}

std::function<double(double, double)>
PotentialFromDensity::
getTruncPotential() const {
    std::function<std::function<std::complex<double>(double, double)>(int, int)>
    wrapper = [=] (int i, int j) {
        return expansion->psi(i, j);
    };
    return utility::realFunction(getTruncFunction(wrapper));
}

std::function<std::array<double, 2>(double, double)>
PotentialFromDensity::
getTruncForce() const {
    std::function<std::function<std::array<std::complex<double>, 2>(double, double)>(int, int)>
    wrapper = [=] (int i, int j) {
        auto debug_output = [i, j, expansion=expansion] (double R, double phi) {
            utility::SimpleTimer timer;
            timer.start();
            // auto time_0 = std::chrono::steady_clock::now();
            auto output = expansion->psi_f(i, j)(R, phi);
            timer.stop();
            utility::debug_print("Force BFE term: " + std::to_string(i) + ", "
                                 + std::to_string(j) + ": " 
                                 + std::to_string(timer.getDuration_us()) + "us", 1);
            // auto time_1 = std::chrono::steady_clock::now();
            // pfd_mtx.lock();
            // std::cout << "Force BFE term: " << i << ", " << j << ", "
            //         << std::chrono::duration_cast<std::chrono::microseconds>
            //                 (time_1 - time_0).count() << "us"
            //         << std::endl;
            // pfd_mtx.unlock();
            return output;
        };
        return debug_output;

        // return expansion->psi_f(i, j); // usually this is all
    };

    // return utility::realFunction(getTruncFunction(wrapper));
    
    return [this, wrapper] (double R, double phi) {
        utility::SimpleTimer timer;
        timer.start();
        auto output = utility::realFunction(getTruncFunction(wrapper))(R, phi);
        timer.stop();
        utility::debug_print("Trunc force: " + std::to_string(timer.getDuration_us())
                             + "us", 1);
        return output;
    }
}

std::vector<std::vector<std::complex<double>>>
PotentialFromDensity::getCoefficients() const {
    return coefficients;
}

double PotentialFromDensity::calculateAbsNorm() const {
    double output = 0;
    for(int i = 0; i <= nl_max[0]; ++i) {
        for(int j = 0; j <= nl_max[1]; ++j) {
            output += pow(std::abs(coefficients[i][j]), 2);
        }
    }
    return sqrt(output);
}

double PotentialFromDensity::
getParticleMass(const std::vector<std::array<double, 2>> positions,
                const std::function<double(double, double)>
                target_density, 
                const std::array<double, 2> target_coords) {
    int N_particles = positions.size();
    std::vector<std::array<double, 3>> density(N_particles);
    for(int i = 0; i < N_particles; ++i) {
        density[i] = {1, positions[i][0], positions[i][1]};
    }
    initFromDensity(density);
    double mass = target_density(target_coords[0], target_coords[1])
                  /trunc_density(target_coords[0], target_coords[1]);
    return mass;
}

double getParticleMass(const std::vector<std::array<double, 2>> positions,
                       const std::function<double(double, double)>
                       target_density, 
                       const std::array<double, 2> target_coords) {
    PotentialFromDensity pfd;
    return pfd.getParticleMass(positions, target_density, target_coords);
}

template void PotentialFromDensity::
         initFromDensity(const std::function<double(double, double)>& density);
template void PotentialFromDensity::
         initFromDensity(const std::vector<std::array<double, 3>>& density);
template std::vector<std::vector<std::complex<double>>>
         PotentialFromDensity::
         getCoefficients(const std::function<double(double, double)>& density) const;
template std::vector<std::vector<std::complex<double>>>
         PotentialFromDensity::
         getCoefficients(const std::vector<std::array<double, 3>>& density) const;
template std::vector<std::vector<std::function<std::complex<double>(double, double)>>>
         PotentialFromDensity::
         getTerms(std::function<std::function<std::complex<double>(double, double)>
                    (int, int)> BFE_member_function) const;
template std::vector<std::vector<std::function<std::array<std::complex<double>, 2>(double, double)>>>
         PotentialFromDensity::
         getTerms(std::function<std::function<std::array<std::complex<double>, 2>(double, double)>
                    (int, int)> BFE_member_function) const;
template std::function<std::complex<double>(double, double)>
        PotentialFromDensity::
        getTruncFunction(std::function<std::function<std::complex<double>(double, double)>
                    (int, int)> BFE_member_function) const;
template std::function<std::array<std::complex<double>, 2>(double, double)>
        PotentialFromDensity::
        getTruncFunction(std::function<std::function<std::array<std::complex<double>, 2>(double, double)>
                    (int, int)> BFE_member_function) const;


}