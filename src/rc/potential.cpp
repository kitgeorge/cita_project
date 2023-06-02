/**
 * \file		shared/inc/potential.cpp
 * \brief		Axisymmetric galactic potential.
 * \author	Rimpei Chiba
 * \date		2020-
 */

#include "potential.hpp"

namespace RC
{
	double Potential::pot(const Polar &pos)
	{
		return pot(pos.r);
	}

	Force Potential::force(const Polar &pos)
	{
		const auto fr = - dPhidr(pos.r);
		Force f;
		f.fR = fr * pos.st;
		f.fz = fr * pos.ct;
		return f;
	}

	double Potential::pot(double r)
	{
		return vc_sun_sq * log(r / R_sun);
	}

	double Potential::dPhidr(double r)
	{
		return vc_sun_sq / r;
	}

	double Potential::d2Phidr2(double r)
	{
		return - vc_sun_sq / (r * r);
	}

	double Potential::rcGivenL(double L)
	{
		// Circular radius rc given L
		// L = rc * vc(rc) = rc * sqrt(rc * dPhi/dr)
		// L^2 = rc^3 * dPhi/dr
		// Solve f(rc) = rc^3 * dPhi/dr - L^2 = 0
		// f(rc) has only one root at rc > 0
		const auto drc    = 1.e-0 * Units::kpc;
		const auto rc_err = 1.e-7 * Units::kpc;
		const auto Lsq = L * L;
		double rc = rc_err;
		while(rc * rc * rc * dPhidr(rc) - Lsq < 0) rc += drc;
		rc = BracketNewtonRaphson(rc - drc, rc, [&](double x, double &dfdx) 
		{
			const auto xsq = x * x;
			dfdx = xsq * (3. * dPhidr(x) + x * d2Phidr2(x));
			return xsq * x * dPhidr(x) - Lsq;
		}, rc_err);
		if(rc < 0) rc = 0.;
		return rc;
	}

	double Potential::rcGivenE(double E)
	{
		// Circular radius rc given E
		// E = 1/2 vc(rc)^2 + Phi(rc)
		//   = 1/2 rc * dPhi/dr(rc) + Phi(rc)
		// Solve f(rc) = 1/2 rc * dPhi/dr(rc) + Phi(rc) - E = 0
		// f is an increasing function
		const auto drc    = 1.e-0 * Units::kpc;
		const auto rc_err = 1.e-7 * Units::kpc;
		double rc = rc_err;
		while(0.5 * rc * dPhidr(rc) + pot(rc) - E < 0) rc += drc;
		rc = BracketNewtonRaphson(rc - drc, rc, [&](double x, double &dfdx) 
		{
			dfdx = 0.5 * (dPhidr(x) + x * d2Phidr2(x)) + dPhidr(x);
			return 0.5 * x * dPhidr(x) + pot(x) - E;
		}, rc_err);
		if(rc < 0) rc = 0.;
		return rc;
	}

	double Potential::rcGivenOmega(double Omega)
	{
		// Circular radius rc given circular frequency Omega
		// Omega(rc)^2 = (1 / rc) * dPhi/dr(rc)
		// Solve f(rc) = Omega^2 - Omega(rc)^2 = 0
		//             = Omega^2 - (1 / rc) * dPhi/dr(rc) = 0
		// f is an increasing function
		const auto drc    = 1.e-0 * Units::kpc;
		const auto rc_err = 1.e-7 * Units::kpc;
		const auto wsq = Omega * Omega;
		double rc = rc_err;
		while(wsq - dPhidr(rc) / rc < 0) rc += drc;
		rc = BracketNewtonRaphson(rc - drc, rc, [&](double x, double &dfdx) 
		{
			dfdx = (- d2Phidr2(x) + dPhidr(x) / x) / x;
			return wsq - dPhidr(x) / x;
		}, rc_err);
		if(rc < 0) rc = 0.;
		return rc;
	}

	double Potential::vcGivenrc(double rc)
	{
		// Circular velocity given rc
		// vc(rc) = sqrt(rc * dPhi/dr)
		return sqrt(rc * dPhidr(rc));
	}

	double Potential::LcGivenrc(double rc)
	{
		// Angular momentum given rc
		// Lc = rc * vc(rc) = rc * sqrt(rc * dPhi/dr)
		return rc * vcGivenrc(rc);
	}

	double Potential::EcGivenrc(double rc)
	{
		// Circular energy given rc
		// Ec = 0.5 * vc(rc)^2 + Phi(rc)
		//    = 0.5 * rc * (dPhi/dr)_rc + Phi(rc)
		return 0.5 * rc * dPhidr(rc) + pot(rc);
	}

	double Potential::EcGivenL(double L)
	{
		// Circular energy given L
		return EcGivenrc(rcGivenL(L));
	}

	double Potential::Omega(double rc)
	{
		// Circular frequency Omega(rc) given rc
		// Omega^2(rc) = (1/R dPhi/dr)_rc
		return sqrt(dPhidr(rc) / rc);
	}

	double Potential::kappa(double rc)
	{
		// Epicycle frequency kappa(rc) given rc
		// kappa^2 = (d2Phieff/dr^2)_rc
		//         = (d2Phi/dr^2)_rc + 3L^2/rc^4
		//         = (d2Phi/dr^2)_rc + 3/rc (dPhi/dr)_rc
		return sqrt(d2Phidr2(rc) + 3. / rc * dPhidr(rc));
	}
	
	double Potential::rcGivenwpN(double Omegap, int Nr, int Npsi, int Nphi)
	{
		// f(rc) = Nr*kappa(rc) + Npsi*Omega(rc) - Nphi*Omegap = 0
		// f is generally a decreasing function
		const auto warn = false;
		const auto rc_limit = 5.e+2 * Units::kpc;
		const auto rc_err   = 1.e-7 * Units::kpc;
		auto       drc      = 1.e-0 * Units::kpc;
		auto rc = rc_err;
		if(Nr * kappa(rc_err) + Npsi * Omega(rc_err) - Nphi * Omegap > 0) // If N•Omega decreases with R
		{
			while(Nr * kappa(rc) + Npsi * Omega(rc) - Nphi * Omegap > 0)
			{
				rc += drc;
				if(rc > rc_limit) break;
			}
		}
		else // For ILRs where N•Omega tends to zero towards R=0
		{
			drc = 0.1 * Units::kpc;
			while(Nr * kappa(rc) + Npsi * Omega(rc) - Nphi * Omegap < 0)
			{
				rc += drc;
				if(rc > rc_limit) break;
			}
			if(Nr == -1)
			{
				rc = rc_err;
				std::ofstream fout("test.dat");
				while(Nr * kappa(rc) + Npsi * Omega(rc) - Nphi * Omegap < 0)
				{
					fout << rc << " " << Nr * kappa(rc) + Npsi * Omega(rc) - Nphi * Omegap << std::endl;
					rc += drc;
					if(rc > rc_limit) break;
				}
				fout.close();
			}
		}
		if(rc > rc_limit)
		{
			if(warn) std::cerr << "Warning [shared/src/potential.cpp/rcGivenwpN]: rc > "
				<< rc_limit << ". (Omegap = " << Omegap * Units::Gyr << " [/Gyr], Nr = " 
				<< Nr << ", Npsi = " << Npsi << ", Nphi = " << Nphi << ")" << std::endl;
		}
		else if(rc < 0)
		{
			if(warn) std::cerr << "Warning [shared/src/potential.cpp/rcGivenwpN]: rc > 0. "
				"(Omegap = " << Omegap * Units::Gyr << " [/Gyr], Nr = " << Nr 
				<< ", Npsi = " << Npsi << ", Nphi = " << Nphi << ")" << std::endl;
		}
		else
		{
			rc = bisection(rc - drc, rc, [&](double x) 
			{
				return Nr * kappa(x) + Npsi * Omega(x) - Nphi * Omegap;
			}, rc_err);
		}
		return rc;
	}
}
