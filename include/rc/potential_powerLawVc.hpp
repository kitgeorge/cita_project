/**
 * \file		shared/inc/potential_powerLawVc.hpp
 * \brief		Potentials for power-law rotation curves.
 * \author	Rimpei Chiba
 * \date		2020-
 */

#ifndef _POWERLAWVC_H_
#define _POWERLAWVC_H_

#include <stdlib.h>
#include <cstdio>
#include <math.h>

#include "galaxy.hpp"
#include "potential.hpp"
#include "utils.hpp"

namespace RC
{
	class PowerLawVc: public Potential
	{
		private:
			const double beta, two_beta;
		
		public:
			double pot(const Polar &pos);
			Force force(const Polar &pos);
			double pot(double r);
			double dPhidr(double r);
			double d2Phidr2(double r);
			double rcGivenL(double L);
			double rcGivenE(double E);
			double rcGivenOmega(double Omega);
			double vcGivenrc(double rc);
			double LcGivenrc(double rc);
			double EcGivenrc(double rc);
			double Omega(double rc);
			double kappa(double rc);
			double rcGivenwpN(double Omegap, int Nr, int Npsi, int Nphi);

			PowerLawVc(double beta);
	};
}

#endif
