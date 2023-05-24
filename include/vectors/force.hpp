#pragma once

#include <array>

namespace vectors {

struct Force {

    Force()=default;
    Force(std::array<double, 3> components);
    void set_components(std::array<double, 3> components);

    // f_z not currently in use (2d simulations)

    double f_R;
    double f_phi;
    double f_z;

    void multiply(double scalar);
    void add(Force f);



    Force operator + (const Force& f);
    void operator += (const Force& f);
    void operator -= (const Force& f);
    void operator *= (const double& scalar);
    void operator /= (const double& scalar);


};
}