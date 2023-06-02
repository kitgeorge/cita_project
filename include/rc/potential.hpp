/**
 * \file		shared/inc/potential.hpp
 * \brief		Axisymmetric galactic potential.
 * \author	Rimpei Chiba
 * \date		2020-
 */

#ifndef _RC_POTENTIAL_H_
#define _RC_POTENTIAL_H_

#include <stdlib.h>
#include <math.h>

#include "galaxy.hpp"
#include "force.hpp"
#include "coordinate.hpp"
#include "NewtonRaphson.hpp"

namespace RC
{
  class Potential
  {
    public:
      virtual double pot(const Polar &pos);
      virtual Force force(const Polar &pos);

      // Functions only for spherical potentials.
      // For non-spherical potentials,
      // these are evaluated at the galactic mid-plane (z=0).
      virtual double pot(double r);
      virtual double dPhidr(double r);
      virtual double d2Phidr2(double r);
      virtual double rcGivenL(double L);
      virtual double rcGivenE(double E);
      virtual double rcGivenOmega(double Omega);
      virtual double vcGivenrc(double rc);
      virtual double LcGivenrc(double rc);
      virtual double EcGivenrc(double rc);
      virtual double EcGivenL(double L);
      virtual double Omega(double rc);
      virtual double kappa(double rc);
      virtual double rcGivenwpN(double Omegap, int Nr, int Npsi, int Nphi);
  };
}

#endif
