#include "potential_funcs.hpp"
#include <iostream>
#include <mutex>

std::mutex pf_mtx;

namespace potential {

double PotentialFuncs::potential(double R, double phi, double t) const {
    return utility::addFunctions(potentials, R, phi, t);
}

std::array<double, 2> 
PotentialFuncs::polar_force(double R, double phi, double t) const {
    return utility::addFunctions(polar_forces, R, phi, t);
}

std::array<double, 2> 
PotentialFuncs::cartesian_force(double R, double phi, double t) const {
    return utility::addFunctions(cartesian_forces, R, phi, t);
}

std::function<std::array<double, 2>(double, double, double)>
PotentialFuncs::
getCartesianForce(std::function<std::array<double, 2>(double, double, double)>
                  polar_force) {
    return [polar_force](double x, double y, double t) {
        std::array<std::array<double, 2>, 2>
        coords = {{ {{x, y}}, {{0, 0}} }};
        vrs::Coords2d pos(coords, 0);
        
        pf_mtx.lock();
        std::cout << "getCartesianForce, " << coords[0][0] << ", "
                  << coords[0][1] << ", " << pos.polar[0][0] << ", "
                  << pos.polar[0][1] << std::endl;
        assert(std::isfinite(coords[0][0]));
        assert(std::isfinite(coords[0][1]));
        assert(std::isfinite(pos.polar[0][0]));
        assert(std::isfinite(pos.polar[0][1]));
        pf_mtx.unlock();

        std::array<double, 2> f = {polar_force(pos.polar[0][0], 
                                                pos.polar[0][1], t)[0],
                                    polar_force(pos.polar[0][0], 
                                                pos.polar[0][1], t)[1]};
        f = vrs::getCartesianVector2d(pos.polar[0], f);
        return f;
    };
}

PotentialFuncs::
PotentialFuncs(std::function<double(double, double, double)> potential,
               std::function<std::array<double, 2>(double, double, double)>
               polar_force):
               potentials({potential}), polar_forces({polar_force}) ,
               cartesian_forces({getCartesianForce(polar_force)}) {}        

PotentialFuncs::PotentialFuncs(const basis_functions::PotentialFromDensity& p):
    potentials({[p] (double R, double phi, double t) {
        assert(std::isfinite(R));
        return p.trunc_potential(R, phi);
    }}),
    polar_forces({[p] (double R, double phi, double t) {
        assert(std::isfinite(R));
        return p.trunc_force(R, phi);
    }}),
    cartesian_forces({getCartesianForce(polar_forces[0])}) {}

PotentialFuncs::PotentialFuncs(const PotentialFuncs& old) {
    potentials = old.potentials;
    polar_forces = old.polar_forces;
    cartesian_forces = old.cartesian_forces;
}

PotentialFuncs::PotentialFuncs(const std::vector<PotentialFuncs>& p) {
    int N = p.size();
    for(int i = 0; i < N; ++i) {
        (*this) += p[i];
    }
}

void PotentialFuncs::operator +=(const PotentialFuncs& p) {
    int N = p.potentials.size();
    assert(p.polar_forces.size() == N);
    assert(p.cartesian_forces.size() == N);

    int N_0 = potentials.size();
    assert(polar_forces.size() == N_0);
    assert(cartesian_forces.size() == N_0);

    potentials.resize(N_0 + N);
    polar_forces.resize(N_0 + N);
    cartesian_forces.resize(N_0 + N);
    for(int i = 0; i < N; ++i) {
        potentials[N_0 + i] = p.potentials[i];
        polar_forces[N_0 + i] = p.polar_forces[i];
        cartesian_forces[N_0 + i] = p.cartesian_forces[i];
    }
}

PotentialFuncs PotentialFuncs::operator *(const double& factor) const {
    PotentialFuncs output(*this);
    for(int i = 0; i < output.potentials.size(); ++i) {
        output.potentials[i] = utility::
                multiplyFunction(output.potentials[i], factor);
        output.polar_forces[i] = utility::
                multiplyFunction(output.polar_forces[i], factor);
        output.cartesian_forces[i] = utility::
                multiplyFunction(output.cartesian_forces[i], factor);
    }
    return output;
}


void PotentialFuncs::multiply(std::function<double(double, double, double)>
                              envelope) {
    int N = potentials.size();
    assert(polar_forces.size() == N);
    assert(cartesian_forces.size() == N);

    for(int i = 0; i < N; ++i) {
        std::function<double(double, double, double)>
        potential_0 = potentials[i];
        std::function<std::array<double, 2>(double, double, double)>
        polar_force_0 = polar_forces[i];
        std::function<std::array<double, 2>(double, double, double)>
        cartesian_force_0 = cartesian_forces[i];

        potentials[i] = [potential_0, envelope] (double R, double phi, double t) {
            return envelope(R, phi, t)*potential_0(R, phi, t);
        };

        polar_forces[i] = [envelope, polar_force_0] (double R, double phi, double t) {
            std::array<double, 2> f = polar_force_0(R, phi, t);
            for(int i = 0; i < 2; ++i) {
                f[i] = envelope(R, phi, t)*f[i];
            }
            return f;
        };

        cartesian_forces[i] = [envelope, cartesian_force_0](double x, double y, double t) {
            std::array<double, 2> f = cartesian_force_0(x, y, t);
            double R = sqrt(pow(x, 2) + pow(y, 2));
            double phi = atan(y/x) + std::numbers::pi*(x < 0);
            for(int i = 0; i < 2; ++i) {
                f[i] = envelope(R, phi, t)*f[i];
            }
            return f;
        };
    }
    }
} 