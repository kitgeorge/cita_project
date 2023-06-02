/**
 * \file		shared/inc/potential_Hernquist.cpp
 * \brief		Hernquist model.
 * \author	Rimpei Chiba
 * \date		2020-
 */

#include "potential_Hernquist.hpp"

namespace RC
{
	//
	// Hernquist (potential)
	//

	Hernquist::Hernquist(double rs, double M):
		rs(rs),
		GM(Units::G * M)
	{}

	double Hernquist::pot(const Polar &pos)
	{
		return pot(pos.r);
	}

	Force Hernquist::force(const Polar &pos)
	{
		const auto fr = - dPhidr(pos.r);
		Force f;
		f.fR = fr * pos.st;
		f.fz = fr * pos.ct;
		return f;
	}

	double Hernquist::pot(double r)
	{
		return - GM / (r + rs);
	}

	double Hernquist::dPhidr(double r)
	{
		return GM / sq(r + rs);
	}

	double Hernquist::d2Phidr2(double r)
	{
		return - 2. * GM / pow(r + rs, 3);
	}

	double Hernquist::rcGivenE(double E)
	{
		// Circular radius rc given E
		// E = 1/2 vc(rc)^2 + Phi(rc)
		// E = 1/2 rc GM/(rc + rs)^2 - GM/(rc + rs) 
		// Eps = - E / (GM/rs) = - 1/2 x/(x + 1)^2 + 1/(x + 1) ,  x = rc / rs
		// Eps (x + 1)^2 = -1/2 x + (x + 1)
		// Eps (x + 1)^2 - 1/2 x - 1 = 0
		// x^2 + (2 - 1/2Eps) x + 1 - 1/Eps = 0
		// x = (1/4Eps - 1) + sqrt((1/4Eps - 1)^2 + (1/Eps - 1))
		auto Eps = - E / (GM / rs);
		auto a = 1. / (4. * Eps) - 1.;
		return rs * (a + sqrt(a * a + (1. / Eps - 1.)));
	}

	//
	// Hernquist Distribution Function
	//
	
	HernquistDF::HernquistDF(double rs, double M):
		rs(rs),
		vg(sqrt(Units::G * M / rs))
	{}

	double HernquistDF::f(double Eps)
	{
		// Distribution function f normalized to 1 (not M)
		// Numerically confirmed to be 1 upon integration
		if(Eps <= 0 || 1 < Eps) return 0;
		else if(Eps == 1) return 1;
		else
		{
			auto f = (3.*asin(sqrt(Eps)) + sqrt(Eps*(1.-Eps))*(1.-2.*Eps)
				*(8.*Eps*(Eps-1.)-3.)) / (pow(1.-Eps,2.5)*sqrt(2)*pow(TPi*rs*vg,3));
			if(f < 0) f = 0;
			return f;
		}
	}

	double HernquistDF::dndEps(double Eps)
	{
		// Calculate the differential energy distribution
		// dn/dEps = f * g = 1/4pi * f0 * g0
		// ∫_0^1 dEps dn/dEps = 1 confirmed numerically
		// Eps is the normalized binding energy 
		// Eps = - E / (GM/rs)

		if(Eps == 0) return 3.2; // 16/5
		else if(Eps == 1) return 0.;
		else
		{
			// Distribution function f
			// f = M/(sqrt(2)*pow(TPi*rs*vg,3)) * f0
			auto f0 = (3.*asin(sqrt(Eps)) + sqrt(Eps*(1.-Eps))*(1.-2.*Eps)
				*(8.*Eps*(Eps-1.)-3.)) / pow(1.-Eps,2.5);
			if(f0 < 0) f0 = 0.;

			// Density of state g
			// g = 2*sqrt(2)*sqPi*pow(rs*vg,3)/M * g0
			auto g0 = (3.*(4.*Eps*(2.*Eps-1.)+1.)*acos(sqrt(Eps))
				- sqrt(Eps*(1.-Eps))*(4.*Eps-1.)*(2.*Eps+3.)) / (3.*pow(Eps,2.5));
			if(g0 < 0) g0 = 0.;

			return iFPi * f0 * g0;
		}
	}

	double HernquistDF::dfdEps(double Eps)
	{
		//
		// Calculate the gradient of distribution function (normalized to 1)
		// f = C A(Eps)/B(Eps)
		//              dA/dEps * B - A * dB/dEps
		// df/dEps = C --------------------------
		//                         B^2

		auto C = 1. / (sqrt(2) * pow(TPi*rs*vg, 3));
		auto E = 1. - Eps;
		auto A = 3.*asin(sqrt(Eps)) + sqrt(Eps*E)*(1. - 2.*Eps)*(-8.*Eps*E - 3.);
		auto B = pow(E, 2.5);
		auto dAdEps = 64. * pow(Eps*E, 1.5);
		auto dBdEps = - 2.5 * pow(E, 1.5);
		return C * (dAdEps * B - A * dBdEps) / (B * B);
	}
}

