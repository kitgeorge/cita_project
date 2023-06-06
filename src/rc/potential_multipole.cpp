/**
 * \file		shared/src/potential_multipole.cpp
 * \brief		Spherical galactic potential.
 * \author	Rimpei Chiba
 * \date		2022-
 */
/*
#include "potential_multipole.hpp"

namespace RC
{
	Multipole::Multipole(const std::string in)
	{
		readAGAMA(in, Philm, dPhilmdr);
		// Set units
		for(size_t i = 0; i < Philm.size(); ++i)
		{
			Philm[i][0] *= Units::kpc;
			dPhilmdr[i][0] *= Units::kpc;
			for(size_t j = 1; j < Philm[i].size(); ++j)
			{
				// Multiply by TSPi = sqrt(4pi) due to difference in def of spherical harmonics
				Philm[i][j] *= Units::kmskms * TSPi;
				dPhilmdr[i][j] *= Units::kmskms / Units::kpc * TSPi;
			}
		}
		// Set the radius of the left end
		r_min = Philm[0][0];
		log_r_min = log(r_min);
		// Get step size in log radius
		double dlogr{};
		for(size_t i = 1; i < Philm.size(); ++i)
		{
			dlogr += log(Philm[i][0]) - log(Philm[i-1][0]);
		}
		dlogr /= (Philm.size() - 1);
		//
		// Interpolate potential and force using cubic b-spline
		//
		for(size_t j = 1; j < Philm[0].size(); ++j)
		{
			std::vector<double> Philm_vec{}, dPhilmdr_vec{};
			for(size_t i = 0; i < Philm.size(); ++i)
			{
				Philm_vec.push_back(Philm[i][j]);
				dPhilmdr_vec.push_back(dPhilmdr[i][j]);
			}
			boost::math::interpolators::cardinal_cubic_b_spline<double>
				Philm_sp(Philm_vec.data(), Philm_vec.size(), log(r_min), dlogr);
			boost::math::interpolators::cardinal_cubic_b_spline<double>
				dPhilmdr_sp(dPhilmdr_vec.data(), dPhilmdr_vec.size(), log(r_min), dlogr);
			Philm_spline.push_back(Philm_sp);
			dPhilmdr_spline.push_back(dPhilmdr_sp);
		}
	}

	void Multipole::readAGAMA(const std::string &in, 
		std::vector<std::vector<double>> &Philm, 
		std::vector<std::vector<double>> &dPhilmdr)
	{
		std::ifstream data(in);
		if(!data.is_open())
		{
			std::cout << "Warning [shared/src/potential_multipole.cpp]: file \'" 
				<< in << "\' does not exist." << std::endl;
		}
		Philm.clear();
		dPhilmdr.clear();
		std::string line;
		bool readPhilm = false;
		bool readdPhilmdr = false;
		while(std::getline(data, line))
		{
			if(line.find("#Phi") == 0)
			{
				// read two more line2 before reading Phi
				for(int i = 0; i < 2; ++i) std::getline(data, line);
				readPhilm = true;
			}
			if(readPhilm)
			{
				if(!is_blank(line))
				{
					std::stringstream ss(line);
					double value;
					std::vector<double> vec;
					while(ss >> value) vec.push_back(value);
					Philm.push_back(vec);
				}
				else readPhilm = false; // stop reading
			}
			else if(line.find("#dPhi/dr") == 0)
			{
				// read two more lines before reading dPhi/dr
				for(int i = 0; i < 2; ++i) std::getline(data, line);
				readdPhilmdr = true;
			}
			if(readdPhilmdr)
			{
				if(!is_blank(line))
				{
					std::stringstream ss(line);
					double value;
					std::vector<double> vec;
					while(ss >> value) vec.push_back(value);
					dPhilmdr.push_back(vec);
				}
				else
				{
					readdPhilmdr = false; // stop reading
					return;
				}
			}
		}
	}

	double Multipole::pot(const Polar &pos)
	{
		double pot{}, dPhidlogr{};
		for(size_t i = 0; i < Philm_spline.size(); ++i)
		{
			const auto l = (int)sqrt(i);
			const auto m = (int)i - l - l * l;
			//std::cout << "i = " << i << ", l = " << l << ", m = " << m << std::endl;
			const auto Ylm = boost::math::spherical_harmonic(l, m, pos.theta, 0);
			if(pos.r < r_min)
			{
				pot += Philm_spline[i](log_r_min) * Ylm.real();
				dPhidlogr += Philm_spline[i].prime(log_r_min) * Ylm.real();
			}
			else pot += Philm_spline[i](log(pos.r)) * Ylm.real();
		}
		if(pos.r < r_min) pot += dPhidlogr * (pos.r / r_min - 1.);
		// f = f(x0) + df/dx(x0) (x - x0)
		//   = f(x0) + df/dlogx(x0) dlogx/dx(x0) (x - x0)
		//   = f(x0) + df/dlogx(x0) (x / x0 - 1)
		return pot;
	}

	Force Multipole::force(const Polar &pos)
	{
		double fr{}, ft{}, dfrdlogr{}, dftdlogr{}; // f_r, f_θ
		for(size_t i = 0; i < dPhilmdr_spline.size(); ++i)
		{
			const auto l = (int)sqrt(i);
			const auto m = (int)i - l - l * l;
			const auto Ylm  = boost::math::spherical_harmonic(l, m, pos.theta, 0);
			const auto Ylm1 = boost::math::spherical_harmonic(l, m + 1, pos.theta, 0);
			const auto dYlmdtheta = m * pos.cot * Ylm + sqrt((l - m) * (l + m + 1)) * Ylm1;
			if(pos.r < r_min)
			{
				fr += - dPhilmdr_spline[i](log_r_min) * Ylm.real();
				ft += - Philm_spline[i](log_r_min) * dYlmdtheta.real();
				dfrdlogr += - dPhilmdr_spline[i].prime(log_r_min) * Ylm.real();
				dftdlogr += - Philm_spline[i].prime(log_r_min) * dYlmdtheta.real();
			}
			else
			{
				fr += - dPhilmdr_spline[i](log(pos.r)) * Ylm.real();
				ft += - Philm_spline[i](log(pos.r)) * dYlmdtheta.real();
			}
		}
		if(pos.r < r_min)
		{
			fr += dfrdlogr * (pos.r / r_min - 1.);
			ft += dftdlogr * (pos.r / r_min - 1.);
		}
		ft /= pos.r;
		if(pos.z == 0) ft = 0;
		RC::Force f;
		f.fR = fr * pos.st + ft * pos.ct;
		f.fz = fr * pos.ct - ft * pos.st;
		return f;
	}

	double Multipole::pot(double r)
	{
		double pot{}, dPhidlogr{};
		for(size_t i = 0; i < Philm_spline.size(); ++i)
		{
			const auto l = (int)sqrt(i);
			const auto m = (int)i - l - l * l;
			const auto Ylm = boost::math::spherical_harmonic(l, m, Pih, 0);
			if(r < r_min)
			{
				pot += Philm_spline[i](log_r_min) * Ylm.real();
				dPhidlogr += Philm_spline[i].prime(log_r_min) * Ylm.real();
			}
			else pot += Philm_spline[i](log(r)) * Ylm.real();
		}
		if(r < r_min) pot += dPhidlogr * (r / r_min - 1.);
		return pot;
	}

	double Multipole::dPhidr(double r)
	{
		double dPhidr{}, d2Phidrdlogr{};
		for(size_t i = 0; i < dPhilmdr_spline.size(); ++i)
		{
			const auto l = (int)sqrt(i);
			const auto m = (int)i - l - l * l;
			const auto Ylm = boost::math::spherical_harmonic(l, m, Pih, 0);
			if(r < r_min)
			{
				dPhidr += dPhilmdr_spline[i](log_r_min) * Ylm.real();
				d2Phidrdlogr += dPhilmdr_spline[i].prime(log_r_min) * Ylm.real();
			}
			else
			{
				dPhidr += dPhilmdr_spline[i](log(r)) * Ylm.real();
			}
		}
		if(r < r_min) dPhidr += d2Phidrdlogr * (r / r_min - 1.);
		return dPhidr;
	}

	double Multipole::d2Phidr2(double r)
	{
		double d2Phidr2{}, d3Phidr2dlogr{};
		for(size_t i = 0; i < dPhilmdr_spline.size(); ++i)
		{
			const auto l = (int)sqrt(i);
			const auto m = (int)i - l - l * l;
			const auto Ylm = boost::math::spherical_harmonic(l, m, Pih, 0);
			if(r < r_min)
			{
				d2Phidr2 += dPhilmdr_spline[i].prime(log_r_min) * Ylm.real();
				d3Phidr2dlogr += dPhilmdr_spline[i].double_prime(log_r_min) * Ylm.real();
			}
			else
			{
				d2Phidr2 += dPhilmdr_spline[i].prime(log(r)) * Ylm.real();
			}
		}
		if(r < r_min) d2Phidr2 += d3Phidr2dlogr * (r / r_min - 1.);
		d2Phidr2 /= r; // since dPhidr_prime = d2Phidrdlogr = d2Phidr2 / r
		return d2Phidr2;
	}
}
*/