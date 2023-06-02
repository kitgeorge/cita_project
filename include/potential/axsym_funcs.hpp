#pragma once
#include "potential_funcs.hpp"

namespace vrs = vectors;

namespace potential {
// Axsymmetric, time-independent substruct with extra functions
struct AxsymFuncs : public PotentialFuncs {
    std::function<double(double)> potential_R;
    std::function<double(double)> dPhidR;
    std::function<double(double)> d2PhidR2;
    std::function<double(double)> RcGivenL;
    std::function<double(double)> RcGivenE;
    std::function<double(double)> RcGivenOmega;
    std::function<double(double)> vcGivenRc;
    std::function<double(double)> LcGivenRc;
    std::function<double(double)> EcGivenRc;
    std::function<double(double)> EcGivenL;
    std::function<double(double)> Omega;
    std::function<double(double)> kappa;
    std::function<double(double, double, double, double)> RcGivenwpN;

    AxsymFuncs(std::function<double(double, double, double)> potential_,
               std::function<std::array<double, 2>(double, double, double)>
               polar_force_, 
               std::array<std::function<double(double)>, 11> axsym_funcs,
               std::function<double(double, double, double, double)> 
               RcGivenwpN_):
               potential(potential_), polar_force(polar_force_),
               potential_R([potential_] (double R) {return potential_R(R, 0, 0)}),
               dPhidR(axsym_funcs[0]), d2PhidR2(axsym_funcs[1]), 
               RcGivenL(axsym_funcs[2]), RcGivenE(axsym_funcs[3]), 
               RcGivenOmega(axsym_funcs[4]), vcGivenRc(axsym_funcs[5]),
               LcGivenRc(axsym_funs[6]), EcGivenRc(axsym_funcs[7]), 
               EcGivenL(axsym_funcs[8]), Omega(axsym_funcs[9]), 
               kappa(axsym_funcs[10]), RcGivenwpN(RcGivenwpN_) {
    
        cartesian_force = [this](double x, double y, double t) {
            std::array<std::array<double, 2>, 2>
            coords = {{ {{x, y}}, {{0, 0}} }};
            vrs::Coords2d pos(coords, 0);
            std::array<double, 2> f = {polar_force(pos.polar[0][0], 
                                                   pos.polar[0][1], t)[0],
                                       polar_force(pos.polar[0][0], 
                                                   pos.polar[0][1], t)[1]};
            f = vrs::getCartesianVector2d(pos.polar[0], f);
            return f;
        };

    }

    AxsymFuncs(const AxsymFuncs& old):
               potential(old.potential), polar_force(old.polar_force),
               potential_R(old.potential_R), dPhidR(old.dPhidR), 
               d2PhidR2(old.d2PhidR2), RcGivenL(old.RcGivenL), 
               RcGivenE(old.RcGivenE), RcGivenOmega(old.RcGivenOmega),
               vcGivenRc(old.vcGivenRc), LcGivenRc(old.LcGivenRc),
               EcGivenRc(old.EcGivenRc), EcGivenL(old.EcGivenL),
               Omega(old.Omega), kappa(old.kappa), RcGivenwpN(old.RcGivenwpN) {}; 
}

}

