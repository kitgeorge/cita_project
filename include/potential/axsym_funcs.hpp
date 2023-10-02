#pragma once
#include "potential_funcs.hpp"
#include "units.hpp"
#include "NewtonRaphson.hpp"

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
	std::function<double(int, int, double)> resonanceRadius;

	/**
	 * Calculates E under this potential for given coordinates
	 *
	 * @param coords polar coordinates {{R, phi}, {v_R, v_phi}}
	 * @return E value (per unit mass)
	 */
	double EGivenPolar(std::array<std::array<double, 2>, 2> coords);
	/** 
	 * Calculates L under this potential given coordinates
	 *
	 * @param coords polar coordinates {{R, phi}, {v_R, v_phi}}
	 * @return L value (per unit mass)
	 */
	double LGivenPolar(std::array<std::array<double, 2>, 2> coords);

    AxsymFuncs(std::function<double(double, double, double)> potential_,
               std::function<std::array<double, 2>(double, double, double)>
               polar_force_, 
               std::array<std::function<double(double)>, 11> axsym_funcs,
			   std::function<double(int, int, double)> resonanceRadius_):
               PotentialFuncs(potential_, polar_force_), 
               potential_R([potential_] (double R) {return potential_(R, 0, 0);}),
               dPhidR(axsym_funcs[0]), d2PhidR2(axsym_funcs[1]), 
               RcGivenL(axsym_funcs[2]), RcGivenE(axsym_funcs[3]), 
               RcGivenOmega(axsym_funcs[4]), vcGivenRc(axsym_funcs[5]),
               LcGivenRc(axsym_funcs[6]), EcGivenRc(axsym_funcs[7]), 
               EcGivenL(axsym_funcs[8]), Omega(axsym_funcs[9]), 
               kappa(axsym_funcs[10]), resonanceRadius(resonanceRadius_) {
        RcGivenwpN = [this] (double Omegap, int Nr, int Npsi, int Nphi) {
        // Taken from Rimpei's code (potential.cpp)

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
				// std::ofstream fout("test.dat");
				// while(Nr * kappa(rc) + Npsi * Omega(rc) - Nphi * Omegap < 0)
				// {
				// 	fout << rc << " " << Nr * kappa(rc) + Npsi * Omega(rc) - Nphi * Omegap << std::endl;
				// 	rc += drc;
				// 	if(rc > rc_limit) break;
				// }
				// fout.close();
			}
		}
		if(rc > rc_limit)
		{
			// if(warn) std::cerr << "Warning [shared/src/potential.cpp/rcGivenwpN]: rc > "
			// 	<< rc_limit << ". (Omegap = " << Omegap * consts::Gyr << " [/Gyr], Nr = " 
			// 	<< Nr << ", Npsi = " << Npsi << ", Nphi = " << Nphi << ")" << std::endl;
		}
		else if(rc < 0)
		{
			// if(warn) std::cerr << "Warning [shared/src/potential.cpp/rcGivenwpN]: rc > 0. "
			// 	"(Omegap = " << Omegap * consts::Gyr << " [/Gyr], Nr = " << Nr 
			// 	<< ", Npsi = " << Npsi << ", Nphi = " << Nphi << ")" << std::endl;
		}
		else
		{
			rc = bisection(rc - drc, rc, [&](double x) 
			{
				return Nr * kappa(x) + Npsi * Omega(x) - Nphi * Omegap;
			}, rc_err);
		}
		return rc;
	};

    }

    AxsymFuncs(const AxsymFuncs& old):
               PotentialFuncs(old),
               potential_R(old.potential_R), dPhidR(old.dPhidR), 
               d2PhidR2(old.d2PhidR2), RcGivenL(old.RcGivenL), 
               RcGivenE(old.RcGivenE), RcGivenOmega(old.RcGivenOmega),
               vcGivenRc(old.vcGivenRc), LcGivenRc(old.LcGivenRc),
               EcGivenRc(old.EcGivenRc), EcGivenL(old.EcGivenL),
               Omega(old.Omega), kappa(old.kappa), RcGivenwpN(old.RcGivenwpN),
			   resonanceRadius(old.resonanceRadius) {}; 
};

};

