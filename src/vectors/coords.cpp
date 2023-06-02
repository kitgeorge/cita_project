#include "coords.hpp"

namespace vectors {

void Coords2d::updateCartesian() {
    cartesian[0][0] = polar[0][0]*cos(polar[0][1]);
    cartesian[0][1] = polar[0][0]*sin(polar[0][1]);
    cartesian[1] = getCartesianVector2d(polar[0], polar[1]);
}

void Coords2d::updatePolar() {
    polar[0][0] = sqrt(pow(cartesian[0][0], 2) + pow(cartesian[0][1], 2));
    polar[0][1] = atan(cartesian[0][1]/cartesian[0][0]);
    if(cartesian[0][0] < 0) {
        polar[0][1] += std::numbers::pi;
    }
    polar[1] = getPolarVector2d(cartesian[0], cartesian[1]);
}

void Coords2d::setPolar(std::array<std::array<double, 2>, 2> coords) {
    polar = coords;
    updateCartesian();
}

void Coords2d::setCartesian(std::array<std::array<double, 2>, 2> coords) {
    cartesian = coords;
    updatePolar();
}

Coords2d::Coords2d(std::array<std::array<double, 2>, 2> coords, int is_polar) {
    if(is_polar) {
        setPolar(coords);
    }
    else {
        setCartesian(coords);
    }
}

Coords2d::Coords2d(const Coords2d& old) {
    cartesian = old.cartesian;
    polar = old.polar;
}

std::vector<double>
getCartesiansFlat(std::vector<Coords2d> data) {
    std::vector<double> output(data.size()*4);
    for(int i = 0; i < data.size(); ++i) {
        for(int j = 0; j < 2; ++j) {
            for(int k = 0; k < 2; ++k) {
                output[i*4 + j*2 + k] = data[i].cartesian[j][k];
            }
        }
    }
    return output;
}

std::vector<double>
getPolarsFlat(std::vector<Coords2d> data) {
    std::vector<double> output(data.size()*4);
    for(int i = 0; i < data.size(); ++i) {
        for(int j = 0; j < 2; ++j) {
            for(int k = 0; k < 2; ++k) {
                output[i*4 + j*2 + k] = data[i].polar[j][k];
            }
        }
    }
    return output;
}



std::array<double, 2> getCartesianVector2d(std::array<double, 2> pos,
                                           std::array<double, 2> vector) {
    std::array<double, 2> 
    output = {vector[0]*cos(pos[1]) - vector[1]*sin(pos[1]),
              vector[0]*sin(pos[1]) + vector[1]*cos(pos[1])};
    return output;
}

std::array<double, 2> getPolarVector2d(std::array<double, 2> pos,
                                       std::array<double, 2> vector) {
    double phi = atan(pos[1]/pos[0]);
    if(pos[0] < 0) {
        phi += std::numbers::pi;
    }
    std::array<double, 2>
    output = {vector[0]*cos(phi) + vector[1]*sin(phi),
              vector[1]*cos(phi) - vector[0]*sin(phi)};
    return output;
}

}
