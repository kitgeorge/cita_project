/**
 * \file		shared/inc/spiral.cpp
 * \brief		Functions for spiral model.
 * \author	Rimpei Chiba
 * \date		2021
 */

#include "spiral.hpp"

namespace RC
{
  //
  // Spiral potential (R exp(-R))
  //
  // Phi(x,t) = - A * vc^2 / m * Phia(R) * cos[m*p(φ,R,t)]
  // Phia(R)  = R / R0 * exp(- (R - R0) / Rs)
  // p(φ,R,t) = φ - φ0 + wp*(t-t0) + cot(α)ln(R/R0)
  //
  // fR      = - dPhi/dR
  //         = - A * vc^2 / m * [dPhia/dR*cos(m*p) - Phia*m*dp/dR*sin(m*p)]
  // fφ      = - 1/R dPhi/dφ
  //         = - A * vc^2 / R * Phia(R) * sin(m*p)
  //
  // dPhia/dR = 1/R0*exp(-(R-R0)/Rs) + R/R0*(-1/Rs)*exp(-(R-R0)/Rs)
  //          = 1/R0*(1-R/Rs)*exp(-(R-R0)/Rs)
  // dp/dR    = cot(α)/R
  //
  Spiral::Spiral(int m, double A, double Rs, double pitch, 
    double phi0, double Omegap, double t0, double ts):
    m(m),
    A(A),
    Rs(Rs),
    pitch(pitch),
    phi0(phi0),
    Omegap(Omegap),
    t0(t0),
    ts_sq_2(2 * ts * ts)
  {}

  double Spiral::pot(double R, double phi, double t)
  {
    const auto x = R / R_sun;
    const auto dt = t - t0;
    const auto p = (phi - phi0) - Omegap * dt + log(x) / tan(pitch);
    const auto Asp = A * exp(- dt * dt / ts_sq_2);
    return - Asp * vc_sun_sq / m * x * exp(- (R - R_sun) / Rs) * cos(m * p);
  }

  RC::Force Spiral::force(double R, double phi, double t)
  {
    const auto x = R / R_sun;
    const auto dt = t - t0;
    const auto cot_a = 1. / tan(pitch);
    const auto p = (phi - phi0) - Omegap * dt + log(x) * cot_a;
    const auto Asp = A * exp(- dt * dt / ts_sq_2);
    const auto expR = exp(- (R - R_sun) / Rs);
    const auto Phia = x * expR;
    const auto dPhiadR = (1. - R / Rs) / R_sun * expR;
    const auto dpdR = cot_a / R;

    RC::Force f;
    f.fR   = Asp * vc_sun_sq / m * (dPhiadR * cos(m * p) - Phia * m * dpdR * sin(m * p));
    f.fphi = - Asp * vc_sun_sq / R * Phia * sin(m * p);
    return f;
  }


  //
  // Spiral potential (R^2 exp(-R))
  //
  // Phi(x,t) = - A * vc^2 / m * Phia(R) * cos[m*p(φ,R,t)]
  // Phia(R)  = (R / R0)^2 * exp(- (R - R0) / Rs)
  // p(φ,R,t) = φ - φ0 + wp*(t-t0) + cot(α)ln(R/R0)
  //
  // fR      = - dPhi/dR
  //         = - A * vc^2 / m * [dPhia/dR*cos(m*p) - Phia*m*dp/dR*sin(m*p)]
  // fφ      = - 1/R dPhi/dφ
  //         = - A * vc^2 / R * Phia(R) * sin(m*p)
  //
  // dPhia/dR = 2R/R0^2*exp(-(R-R0)/Rs)+(R/R0)^2*(-1/Rs)*exp(-(R-R0)/Rs)
  //          = R/R0^2*(2-R/Rs)*exp(-(R-R0)/Rs)
  // dp/dR    = cot(α)/R
  //

  double Spiral2::pot(double R, double phi, double t)
  {
    const auto x = R / R_sun;
    const auto dt = t - t0;
    const auto p = (phi - phi0) - Omegap * dt + log(x) / tan(pitch);
    const auto Asp = A * exp(- dt * dt / ts_sq_2);
    return - Asp * vc_sun_sq / m * x * x * exp(- (R - R_sun) / Rs) * cos(m * p);
  }

  RC::Force Spiral2::force(double R, double phi, double t)
  {
    const auto x = R / R_sun;
    const auto dt = t - t0;
    const auto cot_a = 1. / tan(pitch);
    const auto p = (phi - phi0) - Omegap * dt + log(x) * cot_a;
    const auto Asp = A * exp(- dt * dt / ts_sq_2);
    const auto expR = exp(- (R - R_sun) / Rs);
    const auto Phia = x * x * expR;
    const auto dPhiadR = x * (2. - R / Rs) / R_sun * expR;
    const auto dpdR = cot_a / R;

    RC::Force f;
    f.fR   = Asp * vc_sun_sq / m * (dPhiadR * cos(m * p) - Phia * m * dpdR * sin(m * p));
    f.fphi = - Asp * vc_sun_sq / R * Phia * sin(m * p);
    return f;
  }
}
