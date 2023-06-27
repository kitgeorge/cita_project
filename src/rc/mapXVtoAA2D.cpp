/**
 * \file		shared/inc/mapXVtoAA2D.cpp
 * \brief		Map (x,v) to (theta,J)
 * \author	Rimpei Chiba
 * \date		2019-
 */


#include "mapXVtoAA2D.hpp"
namespace RC {

MapXVtoAA2D::MapXVtoAA2D(const potential::AxsymFuncs* Phi):
	warn(true),
	Phi(Phi),
	dr_apsis(1.e-2 * Units::kpc),
	dr_apsis_err(1.e-9 * Units::kpc),
	r_upper_limit(2000. * Units::kpc),
	r_lower_limit(1.e-9 * Units::kpc),
	tau_max(2.5), // default 2. Keep it smaller than 2.5
	L(0.), E(0.),
	r_apo(0.), r_peri(0.),
	Jr(0.), Jpsi(0.),
	Tr(0.), Tpsi(0.),
	wr(0.), wpsi(0.),
	thetar(0.), thetapsi(0.)
{}

void MapXVtoAA2D::clear()
{
	L = 0., E = 0., r_apo = 0., r_peri = 0., Jr = 0., Jpsi = 0.;
	Tr = 0., Tpsi = 0., wr = 0., wpsi = 0., thetar = 0., thetapsi = 0.;
}

bool MapXVtoAA2D::findApsis(double L, double E)
{
	const auto Lsq = L * L;
	const auto rc = Phi->RcGivenL(fabs(L));
	if(rc > r_upper_limit)
	{
		if(warn)
		{
			show("Warning [shared/src/MapXVtoAA2D.cpp/findApsis]: rc > r_upper_limit");
			show("rc",rc);
		}
		return false;
	}
	else if(rc < 0)
	{
		if(warn)
		{
			show("Warning [shared/src/MapXVtoAA2D.cpp/findApsis]: rc < 0");
			show("rc",rc);
		}
		return false;
	}
	else if(rc == 0) // radial orbit
	{
		// Set r_peri to r_lower_limit
		r_peri = r_lower_limit;

		// Find r_apo
		// (1) Rough estimate
		r_apo = 0.;
		while((Phi->potential_R(r_apo) - E) < 0)
		{
			r_apo += dr_apsis;
			if(r_apo > r_upper_limit)
			{
				if(warn) show("Warning [shared/src/MapXVtoAA2D.cpp/findApsis]: r_apo > r_upper_limit (case rc==0)");
				return false;
			}
		}

		// (2) Find accurate root
		r_apo = BracketNewtonRaphson(r_apo - dr_apsis, r_apo,
		[&](double r, double &dfdr)
		{
			dfdr = Phi->dPhidR(r);
			return Phi->potential_R(r) - E;
		}, dr_apsis_err);

		// (3) Ensure pr^2 = E - Phi(r_apo) > 0
		while(E - Phi->potential_R(r_apo) <= 0)
		{
			r_apo -= dr_apsis_err;
			if(r_apo < rc)
			{
				r_apo = rc;
				break;
			}
		}
		if(r_apo > r_lower_limit && r_apo < r_upper_limit) return true;
		else return false;
	}
	else // if (0 < rc <= r_upper_limit)
	{
		// (1) First roughly find the initial value of r_apo and r_peri.
		// At r = r_apsis, pr = 0, so E = L^2/2r^2 + Phi(r)
		// To find r_apo (or r_peri), start from rc and increase (or decrease) r
		// until the energy exceeds the required E.

		// Apocentre
		r_apo = rc;
		while(0.5 * Lsq / (r_apo * r_apo) + (Phi->potential_R(r_apo) - E) < 0)
		{
			r_apo += dr_apsis;
			if(r_apo > r_upper_limit)
			{
				if(warn)
				{
					show("Warning [shared/src/MapXVtoAA2D.cpp/findApsis]: r_apo > r_upper_limit");
					show("E", E * Units::kpc2Gyr2_i);
					show("L", L * Units::kpc2Gyr_i);
					show("rc", rc * Units::kpc_i);
					show("r_apo", r_apo * Units::kpc_i);
					show("E_apo", (0.5 * Lsq / (r_apo * r_apo) + Phi->potential_R(r_apo)) * Units::kpc2Gyr2_i);
				}
				return false;
			}
		}
		if(r_apo <= dr_apsis)
		{
			if(warn)
			{
				show("Warning [shared/src/MapXVtoAA2D.cpp/findApsis]: r_apo <= dr_apsis");
				show("rc", rc * Units::kpc_i);
				show("r_apo", r_apo * Units::kpc_i);
			}
			return false;
		}
		if(0.5 * Lsq / ((r_apo - dr_apsis) * (r_apo - dr_apsis)) + (Phi->potential_R(r_apo - dr_apsis) - E) > 0)
		{
			// Near circular orbit. Pass to epicycle approx.
			r_apo = rc;
			return false;
		}

		// Pericentre
		r_peri = rc;
		while(0.5 * Lsq / (r_peri * r_peri) + (Phi->potential_R(r_peri) - E) < 0)
		{
			r_peri -= dr_apsis;
			if(r_peri < r_lower_limit)
			{
				r_peri = r_lower_limit;
				break;
			}
		}
		if(0.5 * Lsq / ((r_peri + dr_apsis) * (r_peri + dr_apsis)) + (Phi->potential_R(r_peri + dr_apsis) - E) > 0)
		{
			// Near circular orbit. Pass to epicycle approx.
			r_peri = rc;
			return false;
		}

		// (2) Find the accurate root of the equation
		// f(r) = Phi(r) + L^2/(2 r^2) - E = 0
		// using Newton-Raphson method.
		r_apo = BracketNewtonRaphson(r_apo - dr_apsis, r_apo,
		[&](double r, double &dfdr)
		{
			auto rsq = r * r;
			dfdr = Phi->dPhidR(r) - Lsq / (rsq * r);
			return Phi->potential_R(r) + 0.5 * Lsq / rsq - E;
		}, dr_apsis_err);
		if(r_apo > r_upper_limit)
		{
			if(warn) show("Warning [shared/src/MapXVtoAA2D.cpp/findApsis]: r_apo > r_upper_limit (after NewtonRaphson)");
			return false;
		}
		if(r_apo < r_lower_limit)
		{
			if(warn) show("Warning [shared/src/MapXVtoAA2D.cpp/findApsis]: r_apo < r_lower_limit (after NewtonRaphson)");
			return false;
		}
		if(r_peri > r_lower_limit)
		{
			r_peri = BracketNewtonRaphson(r_peri, r_peri + dr_apsis,
			[&](double r, double &dfdr)
			{
				auto rsq = r * r;
				dfdr = Phi->dPhidR(r) - Lsq / (rsq * r);
				return Phi->potential_R(r) + 0.5 * Lsq / rsq - E;
			}, dr_apsis_err);
			if(r_peri > r_upper_limit)
			{
				if(warn) show("Warning [shared/src/MapXVtoAA2D.cpp/findApsis]: r_peri > r_upper_limit (after NewtonRaphson)");
				return false;
			}
			if(r_peri < r_lower_limit)
			{
				if(warn) show("Warning [shared/src/MapXVtoAA2D.cpp/findApsis]: r_peri < r_lower_limit (after NewtonRaphson)");
				return false;
			}
		}

		// (3) Adjust r_peri and r_apo such that r_peri_true < r_peri, r_apo < r_apo_true
		// to ensure pr(r) is real (pr^2 > 0) and non-zero within [r_peri,r_apo].
		while((2. * (E - Phi->potential_R(r_apo)) - Lsq / (r_apo * r_apo)) <= 0)
		{
			r_apo -= dr_apsis_err;
			if(r_apo < rc)
			{
				r_apo = rc;
				break;
			}
		}
		while((2. * (E - Phi->potential_R(r_peri)) - Lsq / (r_peri * r_peri)) <= 0)
		{
			r_peri += dr_apsis_err;
			if(r_peri > rc)
			{
				r_peri = rc;
				break;
			}
		}
		//r_apo -= dr_apsis_err;
		//r_peri += dr_apsis_err;
		assert(2*(E - Phi->potential_R(r_apo))  - Lsq/(r_apo*r_apo)   > 0);
		assert(2*(E - Phi->potential_R(r_peri)) - Lsq/(r_peri*r_peri) > 0);
		return true;
	}
}

void MapXVtoAA2D::epicycle(double L)
{  
	// Compute orbital frequencies using epicycle theory.
	const auto rc = Phi->RcGivenL(fabs(L));
	r_apo  = rc;
	r_peri = rc;
	wr     = Phi->kappa(rc);
	wpsi   = Phi->Omega(rc);
	const auto Ec = Phi->EcGivenL(L);
	Jr     = E > Ec ? (E - Ec) / wr : 0;
	Jpsi   = L;
	Tr     = TPi / wr;
	Tpsi   = TPi / wpsi;
	if(L < 0)
	{
		wpsi = - wpsi;
		wr = - wr;
	}
}

bool MapXVtoAA2D::mapLEtoJ(double L_, double E_, int N_tau)
{
	if(isnan(L_) || isnan(E_))
	{
		if(warn)
		{
			show("Warning [shared/src/mapXVtoAA2D.cpp/mapLEtoJ]: L or E is nan.");
			show("L", L_);
			show("E", E_);
		}
		return false;
	}
	clear();
	L = L_, E = E_, Jpsi = L_;
	const auto Lsq  = L * L;
	const auto Emin = Phi->EcGivenL(L);
	if(E <= Emin)
	{
		if(warn)
		{
			show("Warning [shared/src/mapXVtoAA2D.cpp/mapLEtoJ]: E <= Emin. Switch to epicycle theory.");
			show("L   ", L * Units::kpc2Gyr_i);
			show("E   ", E * Units::kpc2Gyr2_i);
			show("Emin", Emin * Units::kpc2Gyr2_i);
		}
		epicycle(L);
		return true;
	}
	else
	{
		// Compute Jr and orbital frequencies by numerical integration.
		// First calculate the apsis (r_apo and r_peri).
		bool success = findApsis(L, E);
		if(!success)
		{
			// Sometimes happens when E > Emin but E - Emin too small
			// or if R_apsis is out of range.
			if(warn)
			{
				show("Warning [shared/src/mapXVtoAA2D.cpp/mapLEtoJ]: findApsis did not converge.");
				show("L   ", L * Units::kpc2Gyr_i);
				show("E   ", E * Units::kpc2Gyr2_i);
				show("Emin", Emin * Units::kpc2Gyr2_i);
			}
			if(fabs((E / Emin - 1)) < 0.0001)
			{
				epicycle(L);
				return true;
			}
			else return false;
		}
		// If success but r_apo and r_peri are too close, switch to epicycle theory
		if(0.5 * fabs(r_apo - r_peri) < dr_apsis_err)
		{
			epicycle(L);
			return true;
		}
		else if(r_apo < r_peri)
		{
			if(warn)
			{
				show("Warning [shared/src/mapXVtoAA2D.cpp/mapLEtoJ]: r_apo < r_peri.");
				show("r_apo ", r_apo * Units::kpc_i);
				show("r_peri", r_peri * Units::kpc_i);
			}
			return false;
		}
		else
		{
			//
			// Calculate Jr, Tr, wr, wpsi by integrating r from 
			// r_peri to r_apo using tanh-sinh quadrature.
			//
			// [tanh-sinh quadrature]
			//     r+             1
			// I = ∫ f(r) dr = Δr ∫ f(r_mid + Δr * x) dx
			//     r-            -1
			//
			// r+ = apoapsis
			// r- = periapsis
			// r_mid = (r+ + r-) / 2
			// Δr    = (r+ - r-) / 2
			//
			//        ∞
			// I = Δr ∫ f(r_mid + Δr * x(τ)) x'(τ) dτ
			//       -∞
			//           N
			//   ≈ Δr Δτ ∑ f(r_mid + Δr * x(Δτ * k)) x'(Δτ * k)
			//          k=-N
			//
			// x(τ) = tanh(π/2 sinh(τ))
			//
			//            π/2 cosh(τ)
			// x'(τ) = ------------------
			//         cosh^2(π/2 sinh(τ))
			//

			assert(r_apo > r_peri);
			// Initialize
			Jr = 0., Tr = 0., Tpsi = 0., wr = 0., wpsi = 0.;
			const auto delta_tau = tau_max / N_tau;
			const auto r_mid   = 0.5 * (r_apo + r_peri);
			const auto delta_r = 0.5 * (r_apo - r_peri);
			const auto delta_r_tau = delta_r * delta_tau;
			bool success = true;
			for(int k = - N_tau; k <= N_tau; ++k)
			{
				const auto tau  = delta_tau * k;
				const auto a    = Pih * sinh(tau);
				const auto g    = tanh(a);
				const auto dgdt = Pih * cosh(tau) / RC::sq(cosh(a));
				const auto r    = r_mid + delta_r * g;
				const auto pr2  = 2. * (E - Phi->potential_R(r)) - Lsq / (r * r);
				if(pr2 > 0)
				{
					const auto pr = sqrt(pr2);
					Jr   += dgdt * pr;
					Tr   += dgdt / pr;
					wpsi += dgdt / (pr * r * r);
				}
				else if(warn && (pr2 < -1.e-9))
				{
					show("Warning [shared/src/mapXVtoAA2D.cpp/mapLEtoJ]: detected pr < 0. Switch to epicycle theory.");
					std::cout << std::setprecision(18)
					<<   "r_peri = " << r_peri * Units::kpc_i
					<< "\nr_apo  = " << r_apo * Units::kpc_i
					<< "\nr      = " << r * Units::kpc_i
					<< "\npr2    = " << pr2 * Units::kpc2Gyr2_i
					<< std::setprecision(6) << std::endl;
					success = false;
					break;
				}
			}
			if(success)
			{
				Jr   *= iPi * delta_r_tau;
				Tr   *= 2. * delta_r_tau;
				wpsi *= 2. * L / Tr * delta_r_tau;
				wr    = L > 0 ? TPi / Tr : - TPi / Tr;
				Tpsi  = TPi / fabs(wpsi);
				return true;
			}
			else
			{
				epicycle(L);
				return true;
			}
		}
	}
}

bool MapXVtoAA2D::mapXVtoAA(double r, double psi, double vr, double vpsi, int N_tau)
{
	clear();
	L = r * vpsi; // could be either pos or neg
	E = Phi->potential_R(r) + 0.5 * (vr * vr + vpsi * vpsi);
	const auto result = mapLEtoJ(L, E, N_tau);
	if(!result)
	{
		if(warn)
		{
			show("Warning [shared/src/mapXVtoAA2D.cpp/mapXVtoAA]: mapLEtoJ failed.");
			show("r", r);
			show("psi", psi);
			show("vr", vr * Units::kms_i);
			show("vpsi", vpsi * Units::kms_i);
		}
		return false;
	}
	else
	{
		// Compute angle variables by integrating the following equations 
		// using tanh-sinh quadrature (to treat the divergence at the boundary).
		//
		//       r                   r
		// θ_r = ∫ dr/pr, θ_φ = φ + ∫ dr/pr (Ω - L/r^2)
		//       r-                  r-
		//
		if(r - r_peri < dr_apsis_err)
		{
			thetar = 0.;
			thetapsi = psi;
		}
		else if(r_apo - r < dr_apsis_err)
		{
			thetar = Pi;
			thetapsi = psi;
		}
		else
		{
			// Initialize
			thetar = 0.;
			double dpsi = 0.;
			const auto delta_tau = tau_max / N_tau;
			const auto r_mid = 0.5 * (r + r_peri);
			const auto delta_r = 0.5 * (r - r_peri);
			const auto delta_r_tau = delta_r * delta_tau;
			const auto Lsq = L * L;
			bool success = true;
			for(int k = - N_tau; k <= N_tau; ++k)
			{
				const auto tau = delta_tau * k;
				const auto a = Pih * sinh(tau);
				const auto g = tanh(a);
				const auto dgdt = Pih * cosh(tau) / (cosh(a) * cosh(a));
				const auto rp = r_mid + delta_r * g;
				const auto rpsqi = 1. / (rp * rp);
				const auto pr2 = 2. * (E - Phi->potential_R(rp)) - Lsq * rpsqi;
				if(pr2 > 0)
				{
					const auto pr = sqrt(pr2);
					thetar += dgdt / pr;
					dpsi += (L * rpsqi - wpsi) * dgdt / pr;
				}
				else
				{
					if(warn)
					{
						show("Warning [shared/src/mapXVtoAA2D.cpp/mapXVtoAA2D]: detected pr2<0. switch to epicycle theory.");
						std::cout << std::setprecision(18)
						<< "r_peri = " << r_peri * Units::kpc_i
						<< "\nr      = " << r * Units::kpc_i
						<< "\nrp     = " << rp * Units::kpc_i
						<< "\npr2    = " << pr2 * Units::kmskms_i
						<< "\nL      = " << L * Units::kpckms_i
						<< "\nReduce tau_max." 
						<< std::setprecision(6) << std::endl;
					}
					success = false;
					break;
				}
			}
			if(success)
			{
				thetar *= wr * delta_r_tau;
				dpsi *= delta_r_tau;
				if(vr * vpsi >= 0)
				{
					thetapsi = psi - dpsi;
				}
				else
				{
					thetar = - thetar;
					thetapsi = psi + dpsi;
				}
			}
			else
			{
				// Epicycle theory.
				const auto rc = Phi->RcGivenL(fabs(L));
				const auto k = Phi->kappa(rc);
				const auto w = Phi->Omega(rc);
				const auto ra = sqrt(2. * Jr / k);
				thetar = atan2(vr, (- k * (r - rc)));
				thetapsi = PItoPI(psi - 2. * w * ra / (k * rc) * sin(thetar) 
					+ 0.5 * Jr / fabs(L) * sin(2. * thetar));
			}
		}		
		return true;
	}
}

double MapXVtoAA2D::W_k(int kr, int kp, const std::function<double (double)> &Philm, int N_tau)
{
	//
	// Calculate W by integrating r from 
	// r_peri to r_apo using tanh-sinh quadrature.
	//
	// [tanh-sinh quadrature]
	//     r+             1
	// I = ∫ f(r) dr = Δr ∫ f(r_mid + Δr * x) dx
	//     r-            -1
	//
	// r = r_mid + Δr * x
	// dr/dx = Δr  ->  dr = Δr dx
	//
	// r+ = apoapsis
	// r- = periapsis
	// r_mid = (r+ + r-) / 2
	// Δr    = (r+ - r-) / 2
	//
	//        ∞
	// I = Δr ∫ f(r_mid + Δr * x(τ)) x'(τ) dτ
	//       -∞
	//           N
	//   ≈ Δr Δτ ∑ f(r_mid + Δr * x(Δτ * k)) x'(Δτ * k)
	//          k=-N
	//
	// x(τ) = tanh(π/2 sinh(τ))
	//
	//            π/2 cosh(τ)
	// x'(τ) = ------------------
	//         cosh^2(π/2 sinh(τ))
	//
	assert(r_apo >= r_peri);
	if(Jr == 0)
	{
		if(kr == 0) return Philm(0.5 * (r_apo + r_peri));
		else return 0;
	}
	else
	{
		const auto delta_tau = tau_max / N_tau;
		const auto r_mid   = 0.5 * (r_apo + r_peri);
		const auto delta_r = 0.5 * (r_apo - r_peri);
		const auto delta_r_tau = delta_r * delta_tau;
		const auto Lsq = L * L;
		double W{}, thetar{}, thetapsi{};
		for(int k = - N_tau; k <= N_tau; ++k)
		{
			const auto tau  = delta_tau * k;
			const auto a    = Pih * sinh(tau);
			const auto g    = tanh(a);
			const auto dgdt = Pih * cosh(tau) / RC::sq(cosh(a));
			const auto r    = r_mid + delta_r * g;
			const auto pr2  = 2. * (E - Phi->potential_R(r)) - Lsq / (r * r);
			if(pr2 > 0)
			{
				const auto pr = sqrt(pr2);
				const auto dthetardr = wr / pr;
				thetar   += dthetardr * dgdt * delta_r_tau;
				thetapsi += - (L / (r * r) - wpsi) / pr * dgdt * delta_r_tau;
				W += Philm(r) * dthetardr * dgdt * cos(kr * thetar + kp * thetapsi);
			}
			else
			{
				if(warn)
				{
					show("Warning [shared/src/mapXVtoAA2D.cpp/W_k]: detected pr < 0.");
					std::cout << std::setprecision(18)
					<<   "r_peri = " << r_peri * Units::kpc_i
					<< "\nr_apo  = " << r_apo * Units::kpc_i
					<< "\nr      = " << r * Units::kpc_i
					<< "\npr2    = " << pr2 * Units::kpc2Gyr2_i
					<< "\nL      = " << L * Units::kpc2Gyr_i
					<< "\nJr     = " << Jr * Units::kpc2Gyr_i
					<< std::setprecision(6) << std::endl;
				}
			}
		}
		W *= iPi * delta_r_tau;
		return W;
	}
}

void gridEOverJ(const potential::AxsymFuncs *Phi, int N_tau, int N_Jr, int N_L, 
double dJr, double dL, std::vector<std::vector<double>> &E_J, bool prog)
{
	E_J.resize(N_Jr, std::vector<double>(N_L, 0.));

	// Set the tolerance and criterion of iteration.
	// With dJr_err=1.e-6, dE=1.e-5, the iteration converges in 3-5 run
	const auto dJr_err = 1.e-6 * Units::kpc2Gyr;
	const auto dE      = 1.e-5 * Units::kpc2Gyr2;
	const auto DE      = 1.e+4 * Units::kpc2Gyr2;
	const auto L_min   = 1.e-2 * Units::kpc2Gyr;
	const int  itr_max = 100;

	MapXVtoAA2D map(Phi);
	const auto PROGRESS = static_cast<int>(N_Jr / 100);
	if(!PROGRESS) prog = false;
	for(int k = 0; k < N_Jr; ++k)
	{
		if(prog && (k % PROGRESS == 0))
			std::cout << "\r" << k / PROGRESS << "%" << std::flush;
		const auto Jr = dJr * k;
		for(int l = 0; l < N_L; ++l)
		{
			const auto L  = dL * l;
			const auto Ec = (L > L_min) ? Phi->EcGivenL(L) : Phi->EcGivenL(L_min);
			double E = Ec + dE;
			// Initial rough search.
			while(1)
			{
				auto success = map.mapLEtoJ(L, E, N_tau);
				if(success && (map.Jr < Jr)) E += DE;
				else break;
			}
			// Newton-Raphson method.
			// f(E) = Jr(E) - Jr_target
			// dfdE = [Jr(E+dE) - Jr(E)] / dE  (> 0)
			// E   -= f / dfdE
			int itr{};
			while(1)
			{
				auto success = map.mapLEtoJ(L, E, N_tau);
				const auto Jr_0 = map.Jr;
				if(success)
				{
					if(abs(map.Jr - Jr) > dJr_err)
					{
						map.mapLEtoJ(L, E + dE, N_tau);
						E -= ((map.Jr - Jr) / (map.Jr - Jr_0)) * dE;
						if(E < Ec) E = Ec;
					}
					else break;
					++ itr;
					if(itr > itr_max) break;
				}
				else break;
			}
			if(fabs(map.Jr - Jr) > dJr_err)
			{
				show("Warning [shared/src/mapXVtoAA2D.cpp]: Jr did not converge.");
				show("L        ", L * Units::kpc2Gyr_i);
				show("Jr_target", Jr * Units::kpc2Gyr_i);
				show("Jr_mapped", map.Jr * Units::kpc2Gyr_i);
				show("E        ", E * Units::kpc2Gyr2_i);
			}
			if(isnan(E) || isinf(E)) E = +1./0.;
			E_J[k][l] = (E < Ec) ? Ec : E;
		}
	}
	if(prog) std::cout << "\r100%" << std::endl;
}

double XGivenJ(double Jr, double L, double dJr, double dL, const std::vector<std::vector<double>> &X_J)
{
	// Get X (bilinear interpolation in 2D action space)
	const auto N_Jr = X_J.size();
	const auto N_L = X_J[0].size();
	const auto i = floor(Jr / dJr);
	const auto j = floor(L / dL);
	if((i < 0) || ((int)N_Jr <= i) || (j < 0) || ((int)N_L <= j))
	{
		/*
		std::cout << "Warning [shared/src/mapXVtoAA2D.cpp/XGivenJ]: "
		<< "Actions are out of range."
		<< "\nJr     = " << Jr * Units::kpckms_i
		<< "\nL      = " << L * Units::kpckms_i
		<< "\nJr_max = " << dJr * (N_Jr - 1) * Units::kpckms_i
		<< "\nL_max  = " << dL * (N_L - 1) * Units::kpckms_i
		<< std::endl;
		*/
		return 0.;
	}
	else if((i < (int)N_Jr - 1) && (j < (int)N_L - 1))
	{
		const auto Jr0 = dJr * i;
		const auto aJr = (Jr - Jr0) / dJr;
		const auto L0  = dL * j;
		const auto aL  = (L - L0) / dL;
		const auto X00 = X_J[i  ][j  ]; // left bottom in (Jr,L)
		const auto X10 = X_J[i+1][j  ]; // right bottom
		const auto X01 = X_J[i  ][j+1]; // left top
		const auto X11 = X_J[i+1][j+1]; // right top
		const auto Xb  = X00 + aJr * (X10 - X00); // Interpolate from X00 to X10.
		const auto Xt  = X01 + aJr * (X11 - X01); // Interpolate from X01 to X11.
		const auto X   = Xb + aL * (Xt - Xb); // Interpolate from Xb to Xt.
		return X;
	}
	else if((i == (int)N_Jr - 1) && (j < (int)N_L - 1))
	{
		const auto L0  = dL * j;
		const auto aL  = (L - L0) / dL;
		const auto XN0 = X_J[i][j  ]; // bottom
		const auto XN1 = X_J[i][j+1]; // top
		const auto X   = XN0 + aL * (XN1 - XN0); // Interpolate from XN0 to XN1.
		return X;
	}
	else if((i < (int)N_Jr - 1) && (j == (int)N_L - 1))
	{
		const auto Jr0 = dJr * j;
		const auto aJr = (Jr - Jr0) / dJr;
		const auto X0N = X_J[i  ][j]; // left
		const auto X1N = X_J[i+1][j]; // right
		const auto X   = X0N + aJr * (X1N - X0N); // Interpolate from X0N to X1N.
		return X;
	}
	else if((i == (int)N_Jr - 1) && (j == (int)N_L - 1))
	{
		return X_J[i][j];
	}
	else // Not possible
	{
		return 0.;
	}
}
}

