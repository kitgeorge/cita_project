/**
 * \file		shared/inc/potential_powerLawVc.cpp
 * \brief		Potentials for power-law rotation curves.
 * \author	Rimpei Chiba
 * \date		2020-
 */

#include "potential_powerLawVc.hpp"

//
// vc(R)  = vc_sun * (R / R_sun)^beta
// Phi(R) = vc_sun^2 / (2 * beta) * (R / R_sun)^(2 * beta)
// dPhi/dR(R) = vc_sun^2 / R * (R / R_sun)^(2 * beta)
//

namespace RC
{
	PowerLawVc::PowerLawVc(double beta):
		beta(beta),
		two_beta(beta + beta)
	{
		if(beta >= 1)
		{
			show("[Error] shared/inc/PowerLawVc.cpp: beta is larger than 1");
			exit(EXIT_FAILURE);
		}
		if(beta <= -1)
		{
			show("[Error] shared/inc/PowerLawVc.cpp: beta is smaller than -1");
			exit(EXIT_FAILURE);
		}
	}

	double PowerLawVc::pot(const Polar &pos)
	{
		return pot(pos.r);
	}

	Force PowerLawVc::force(const Polar &pos)
	{
		const auto fr = - dPhidr(pos.r);
		Force f;
		f.fR = fr * pos.st;
		f.fz = fr * pos.ct;
		return f;
	}

	double PowerLawVc::pot(double r)
	{
		if(beta == 0) return vc_sun_sq * log(r / R_sun);
		else          return vc_sun_sq / two_beta * pow(r / R_sun, two_beta);
	}

	double PowerLawVc::dPhidr(double r)
	{
		if(beta == 0) return vc_sun_sq / r;
		else          return vc_sun_sq / r * pow(r / R_sun, two_beta);
	}

	double PowerLawVc::d2Phidr2(double r)
	{
		if(beta == 0) return - vc_sun_sq / (r * r);
		else          return (two_beta - 1.) * vc_sun_sq / (r * r) * pow(r / R_sun, two_beta);
	}

	double PowerLawVc::rcGivenL(double L)
	{
		if(beta == 0) return L / vc_sun;
		else          return R_sun * pow(L / (vc_sun * R_sun), 1. / (beta + 1.));
	}

	double PowerLawVc::rcGivenE(double E)
	{
		if(beta == 0) return R_sun * exp(E / vc_sun_sq - 0.5);
		else          return R_sun * pow(2. * E / (vc_sun_sq * (1. + 1. / beta)), 1. / two_beta);
	}

	double PowerLawVc::rcGivenOmega(double Omega)
	{
		if(beta == 0) return vc_sun / Omega;
		else          return R_sun * pow(Omega_sun / Omega, 1. / (1. - beta));
	}

	double PowerLawVc::vcGivenrc(double rc)
	{
		if(beta == 0) return vc_sun;
		else          return vc_sun * pow(rc / R_sun, beta);
	}

	double PowerLawVc::LcGivenrc(double rc)
	{
		if(beta == 0) return rc * vc_sun;
		else          return rc * vc_sun * pow(rc / R_sun, beta);
	}

	double PowerLawVc::EcGivenrc(double rc)
	{
		if(beta == 0) return vc_sun_sq * (0.5 + log(rc / R_sun));
		else          return (0.5*(1.+1./beta))*vc_sun_sq*pow(rc/R_sun, two_beta);
	}

	double PowerLawVc::Omega(double rc)
	{
		// Circular frequency Omega(Rc) given Rc
		// Omega^2(Rc) = (1/r dPhi/dR)_Rc
		//             = v0^2 / R^2 * (R/R0)^(2b)
		//             = [v0 / R * (R/R0)^b]^2
		if(beta == 0) return vc_sun / rc;
		else          return vc_sun / rc * pow(rc / R_sun, beta);
	}

	double PowerLawVc::kappa(double rc)
	{
		// Epicycle frequency kappa(Rc) given Rc
		// kappa^2 = (d2Phieff/dR^2)_Rc
		//         = (d2Phi/dR^2)_Rc + 3L^2/Rc^4
		//         = (d2Phi/dR^2)_Rc + 3/Rc (dPhi/dR)_Rc
		// Phi     = v0^2 / (2b) * (R/R0)^(2b)
		// dPhi/dR = v0^2 / R0 * (R/R0)^(2b-1) = v0^2 / R * (R/R0)^(2b)
		// d2Phi/dR^2 = (2b-1) * v0^2 / R0^2 * (R/R0)^(2b-2)
		//            = (2b-1) * v0^2 / R^2 * (R/R0)^(2b)
		// kappa^2 = (2b-1)*Omega^2 + 3*Omega^2
		//         = 2(b+1)*Omega^2
		return sqrt(2. * (beta + 1.)) * Omega(rc);
	}

	double PowerLawVc::rcGivenwpN(double Omegap, int Nr, int Npsi, int Nphi)
	{
		// Rc = R0 * (Omega0/Omega)^(1/(1-beta))
		if(beta == 0)
			return vc_sun * (sqrt(2.) * Nr + Npsi) / (Nphi * Omegap);
		else
			return R_sun*pow(Omega_sun/(Nphi*Omegap)*(sqrt(2.*(beta+1.))*Nr+Npsi), 1./(1.-beta));
	}
}

