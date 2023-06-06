/**
 * \file		shared/inc/potential_multipole.hpp
 * \brief		Multipole expansion.
 * \author	Rimpei Chiba
 * \date		2022-
 */
 /*

#ifndef _RC_POTENTIAL_MULTIPOLE_H_
#define _RC_POTENTIAL_MULTIPOLE_H_

#include <stdlib.h>
#include <iostream>
#include <math.h>
#include <boost/math/special_functions/spherical_harmonic.hpp>
#include <boost/math/interpolators/cardinal_cubic_b_spline.hpp>

#include "math_RC.hpp"
#include "import.hpp"
#include "galaxy.hpp"
#include "force.hpp"
#include "coordinate.hpp"
#include "potential.hpp"
#include "NewtonRaphson.hpp"

//
// Phi(r,θ,φ)    = ∑_l=0^lmax ∑_m=-l^l Philm(r) Ylm(θ,φ)
// dPhidr(r,θ,φ) = ∑_l=0^lmax ∑_m=-l^l dPhilmdr(r) Ylm(θ,φ)
//

namespace RC
{
  class Multipole : public Potential
  {
    private:
      double r_min, log_r_min;
      std::vector<std::vector<double>> Philm{}, dPhilmdr{};
      std::vector<boost::math::interpolators::cardinal_cubic_b_spline<double>> 
        Philm_spline{}, dPhilmdr_spline{};
      void readAGAMA(const std::string &in, 
        std::vector<std::vector<double>> &Philm, 
        std::vector<std::vector<double>> &dPhilmdr);

    public:
      double pot(const Polar &pos);
      Force force(const Polar &pos);
      double pot(double r);
      double dPhidr(double r);
      double d2Phidr2(double r);
      Multipole(const std::string in);
  };
}

#endif
*/