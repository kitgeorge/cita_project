#pragma once
#include <array>
#include <vector>
#include <cmath>
#include <numbers>

namespace vectors {

class Coords2d {
    void updatePolar();
    void updateCartesian();
    public:
        std::array<std::array<double, 2>, 2> cartesian;
        std::array<std::array<double, 2>, 2> polar;

        void setPolar(std::array<std::array<double, 2>, 2> coords);
        void setCartesian(std::array<std::array<double, 2>, 2> coords);

        Coords2d()=default;
        Coords2d(std::array<std::array<double, 2>, 2> coords, int is_polar);
        Coords2d(const Coords2d& old);

};

std::vector<double>
getCartesiansFlat(std::vector<Coords2d> data);

std::vector<double>
getPolarsFlat(std::vector<Coords2d> data);




std::array<double, 2> getCartesianVector2d(std::array<double, 2> pos,
                                           std::array<double, 2> vector);

std::array<double, 2> getPolarVector2d(std::array<double, 2> pos,
                                       std::array<double, 2> vector);
}