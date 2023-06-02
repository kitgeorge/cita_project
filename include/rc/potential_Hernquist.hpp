/**
 * \file		shared/inc/potential_Hernquist.hpp
 * \brief		Hernquist model.
 * \author	Rimpei Chiba
 * \date		2020-
 */

#ifndef _HERNQUIST_H_
#define _HERNQUIST_H_

#include <stdlib.h>
#include <cstdio>
#include <math.h>

#include "galaxy.hpp"
#include "utils.hpp"
#include "potential.hpp"
#include "NewtonRaphson.hpp"
#include "math_RC.hpp"

/**
 * Hernquist model
 * rho(r) = M / (2 pi) * rs / (r * (r + rs)^3)
 *        = rho0 * rs^4 / (r * (r + rs)^3)
 *        = rho0 / (r/rs * (1 + r/rs)^3)
 * Phi(r) = - GM / (r + rs)
 * Input parameters are rs and M = 2 pi rho0 rs^3 (total mass)
 */

namespace RC
{
	class Hernquist: public Potential
	{
		private:
			const double rs, GM;

		public:
			double pot(const Polar &pos);
			Force force(const Polar &pos);
			double pot(double r);
			double dPhidr(double r);
			double d2Phidr2(double r);
			double rcGivenE(double E);
			Hernquist(double rs, double M);
	};

	class HernquistDF
	{
		private:
			const double rs, vg;

		public:
			double f(double Eps);      // Distribution function (isotropic)
			double dndEps(double Eps); // Differential energy distribution 
			double dfdEps(double Eps); // Gradient of distribution function
			HernquistDF(double rs, double M);
	};
}

#endif
