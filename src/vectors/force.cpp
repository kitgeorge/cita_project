#include "force.hpp"

namespace vectors {
void Force::multiply(double scalar) {
    f_R = scalar*f_R;
    f_phi = scalar*f_phi;
    f_z = scalar*f_z;
}

void Force::add(Force f) {
    f_R += f.f_R;
    f_phi += f.f_phi;
    f_z += f.f_z;
}

Force::Force(std::array<double, 3> components) {
    set_components(components);
}

void Force::set_components(std::array<double, 3> components) {
    f_R = components[0];
    f_phi = components[1];
    f_z = components[2];
}

Force Force::operator + (const Force& f) {
    std::array<double, 3> components = {f_R + f.f_R, f_phi + f.f_phi, f_z + f.f_z};
    Force output(components);
    return output;
}

void Force::operator += (const Force& f) {
    add(f);
}

void Force::operator -= (const Force& f) {
    Force temp = f;
    temp.multiply(-1);
    add(temp);
}

void Force::operator *= (const double& scalar) {
    multiply(scalar);
}

void Force::operator /= (const double& scalar) {
    multiply(1/scalar);
}

}