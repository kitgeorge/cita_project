/**
 * \file		shared/src/potential_CartGrid.cpp
 * \brief		Cartesian grid.
 * \author	Rimpei Chiba
 * \date		2022-
 */

/*
#include "potential_CartGrid.hpp"
namespace RC
{
	CartGrid::CartGrid(const std::string in):
	Phi_makima(0, boost::math::interpolators::makima(std::vector<double>{1,2,3,4}, std::vector<double>{1,2,3,4})),
	dPhidR_makima(0, boost::math::interpolators::makima(std::vector<double>{1,2,3,4}, std::vector<double>{1,2,3,4})),
	dPhidz_makima(0, boost::math::interpolators::makima(std::vector<double>{1,2,3,4}, std::vector<double>{1,2,3,4}))
	{
		readAGAMA(in, Phi, dPhidR, dPhidz);
		//
		// Interpolate potential and force using modified Akima
		//
		for(size_t j = 0; j < Phi[0].size(); ++j)
		{
			std::vector<double> x = R, y{};
			for(size_t i = 0; i < Phi.size(); ++i) y.push_back(Phi[i][j]);
			Phi_makima.push_back(boost::math::interpolators::makima(std::move(x), std::move(y)));
		}
		for(size_t j = 0; j < dPhidR[0].size(); ++j)
		{
			std::vector<double> x = R, y{};
			for(size_t i = 0; i < dPhidR.size(); ++i) y.push_back(dPhidR[i][j]);
			dPhidR_makima.push_back(boost::math::interpolators::makima(std::move(x), std::move(y)));
		}
		for(size_t j = 0; j < dPhidz[0].size(); ++j)
		{
			std::vector<double> x = R, y{};
			for(size_t i = 0; i < dPhidz.size(); ++i) y.push_back(dPhidz[i][j]);
			dPhidz_makima.push_back(boost::math::interpolators::makima(std::move(x), std::move(y)));
		}
	}

	void CartGrid::readAGAMA(const std::string &in, 
		std::vector<std::vector<double>> &Phi, 
		std::vector<std::vector<double>> &dPhidR, 
		std::vector<std::vector<double>> &dPhidz)
	{
		std::ifstream data(in);
		if(!data.is_open())
		{
			std::cout << "Warning [shared/src/potential_CartGrid.cpp]: file \'" 
				<< in << "\' does not exist." << std::endl;
		}
		Phi.clear();
		dPhidR.clear();
		dPhidz.clear();
		std::string line;
		bool readPhi = false;
		bool readdPhidR = false;
		bool readdPhidz = false;
		while(std::getline(data, line))
		{
			//
			// Read Phi (Phi[iz][iR])
			//
			if(line.find("type=CylSpline") == 0)
			{
				// read seven more lines before reading Phi
				for(int i = 0; i < 8; ++i) std::getline(data, line);
				readPhi = true;
			}
			if(readPhi)
			{
				if(!is_blank(line))
				{
					std::stringstream ss(line);
					double value;
					std::vector<double> vec;
					for(; ss; )
					{
						if(ss >> value) vec.push_back(value);
						else if(!ss.eof())
						{
							ss.clear();
							ss.ignore(1);
						}
					}
					Phi.push_back(vec);
				}
				else readPhi = false; // stop reading
			}
			//
			// Read dPhi/dR
			//
			else if(line.find("#dPhi/dR") == 0)
			{
				// read two more lines before reading dPhi/dR
				for(int i = 0; i < 2; ++i) std::getline(data, line);
				readdPhidR = true;
			}
			if(readdPhidR)
			{
				if(!is_blank(line))
				{
					std::stringstream ss(line);
					double value;
					std::vector<double> vec;
					for(; ss; )
					{
						if(ss >> value) vec.push_back(value);
						else break; // ignore first line "z"
					}
					if(vec.size())
					{
						vec.erase(vec.begin()); // erase first element "R"
						dPhidR.push_back(vec);
					}
				}
				else readdPhidR = false; // stop reading
			}
			//
			// Read dPhi/dz
			//
			else if(line.find("#dPhi/dz") == 0)
			{
				// read two more lines before reading dPhi/dz
				for(int i = 0; i < 2; ++i) std::getline(data, line);
				readdPhidz = true;
			}
			if(readdPhidz)
			{
				if(!is_blank(line))
				{
					std::stringstream ss(line);
					double value;
					std::vector<double> vec;
					for(; ss; )
					{
						if(ss >> value) vec.push_back(value);
						else break; // ignore first line "z"
					}
					if(vec.size())
					{
						vec.erase(vec.begin()); // erase first element "R"
						dPhidz.push_back(vec);
					}
				}
				else
				{
					readdPhidz = false; // stop reading
					return;
				}
			}
		}
		// Set z
		z.clear();
		z = Phi[0];
		z_max = z[z.size()-1];
		Phi.erase(Phi.begin());
		// Set R
		R.clear();
		for(size_t i = 0; i < Phi.size(); ++i)
		{
			R.push_back(Phi[i][0]);
			Phi[i].erase(Phi[i].begin());
		}
		R_max = R[R.size()-1];
		//
		// Set units
		//
		for(size_t iz = 0; iz < Phi.size(); ++iz)
		{
			for(size_t iR = 0; iR < Phi[iz].size(); ++iR)
			{
				Phi[iz][iR]    *= Units::kmskms;
				dPhidR[iz][iR] *= Units::kmskms / Units::kpc;
				dPhidz[iz][iR] *= Units::kmskms / Units::kpc;
			}
		}
		for(auto& a : R) a *= Units::kpc;
		for(auto& a : z) a *= Units::kpc;
	}

	double CartGrid::pot(const Polar &pos)
	{
		if((abs(pos.z) < z_max) && (pos.R < R_max))
		{
			size_t iz = 0;
			for(; iz < z.size(); ++iz) if(z[iz] > abs(pos.z)) break;
			std::vector<double> z_{}, Phi_{};
			for(int i = iz - 2; i < (int)iz + 2; ++i)
			{
				if(iz + 2 > Phi_makima.size()) // if z is in the last interval
				{
					z_.push_back(z[i-1]);
					Phi_.push_back(Phi_makima[abs(i-1)](pos.R));
				}
				else
				{
					z_.push_back(i < 0 ? - z[abs(i)] : z[i]);
					Phi_.push_back(Phi_makima[abs(i)](pos.R));
				}
			}
			boost::math::interpolators::makima z_spline(std::move(z_), std::move(Phi_));
			return z_spline(abs(pos.z));
		}
		else
		{
			// Extrapolate from the boundary of the CartGrid assuming a 1/r decay
			double Rb{}, zb{}; // (R,z) at the boundary of the CartGrid
			// If slope in (R,z) is smaller than z_max / R_max
			if((abs(pos.z) * R_max) < (pos.R * z_max))
			{
				Rb = R_max;
				zb = R_max / pos.R * abs(pos.z);
			}
			else
			{
				Rb = z_max / abs(pos.z) * pos.R;
				zb = z_max;
			}
			const auto r_ratio = sqrt((RC::sq(Rb) + RC::sq(zb)) / (RC::sq(pos.R) + RC::sq(pos.z)));

			size_t iz = 0;
			for(; iz < z.size(); ++ iz) if(z[iz] > zb) break;
			std::vector<double> z_{}, Phi_{};
			for(int i = iz - 2; i < (int)iz + 2; ++i)
			{
				if(iz + 2 > Phi_makima.size()) // if z is in the last interval
				{
					z_.push_back(z[i-1]);
					Phi_.push_back(Phi_makima[abs(i-1)](Rb));
				}
				else
				{
					z_.push_back(i < 0 ? - z[abs(i)] : z[i]);
					Phi_.push_back(Phi_makima[abs(i)](Rb));
				}
			}
			boost::math::interpolators::makima z_spline(std::move(z_), std::move(Phi_));
			return z_spline(zb) * r_ratio;
		}
	}

	Force CartGrid::force(const Polar &pos)
	{
		Force f{};
		if((abs(pos.z) < z_max) && (pos.R < R_max))
		{
			size_t iz = 0;
			for(; iz < z.size(); ++iz) if(z[iz] >= abs(pos.z)) break;
			std::vector<double> zR{}, dPhidR_{}, dPhidz_{};
			for(int i = iz - 2; i < (int)iz + 2; ++i)
			{
				if(iz + 2 > dPhidR_makima.size())
				{
					zR.push_back(z[i-1]);
					dPhidR_.push_back(dPhidR_makima[abs(i-1)](pos.R));
					dPhidz_.push_back(dPhidz_makima[abs(i-1)](pos.R));
				}
				else
				{
					zR.push_back(i < 0 ? - z[abs(i)] : z[i]);
					dPhidR_.push_back(dPhidR_makima[abs(i)](pos.R));
					dPhidz_.push_back(dPhidz_makima[abs(i)](pos.R));
				}
			}
			std::vector<double> zz = zR;
			boost::math::interpolators::makima dPhidR_spline(std::move(zR), std::move(dPhidR_));
			boost::math::interpolators::makima dPhidz_spline(std::move(zz), std::move(dPhidz_));
			f.fR = - dPhidR_spline(abs(pos.z));
			f.fz = pos.z < 0 ? dPhidz_spline(abs(pos.z)) : - dPhidz_spline(abs(pos.z));
		}
		else
		{
			// Extrapolate from the boundary of the CartGrid assuming a 1/r^2 decay
			double Rb{}, zb{}; // (R,z) at the boundary of the CartGrid
			// If slope in (R,z) is smaller than z_max / R_max
			if((abs(pos.z) * R_max) < (pos.R * z_max))
			{
				Rb = R_max;
				zb = R_max / pos.R * abs(pos.z);
			}
			else
			{
				Rb = z_max / abs(pos.z) * pos.R;
				zb = z_max;
			}
			const auto r_sq_ratio = (RC::sq(Rb) + RC::sq(zb)) / (RC::sq(pos.R) + RC::sq(pos.z));

			size_t iz = 0;
			if(zb == z_max) iz = z.size() - 1;
			else for(; iz < z.size(); ++ iz) if(z[iz] >= zb) break;
			std::vector<double> zR{}, dPhidR_{}, dPhidz_{};
			for(int i = iz - 2; i < (int)iz + 2; ++i)
			{
				if(iz + 2 > dPhidR_makima.size())
				{
					zR.push_back(z[i-1]);
					dPhidR_.push_back(dPhidR_makima[abs(i-1)](Rb));
					dPhidz_.push_back(dPhidz_makima[abs(i-1)](Rb));
				}
				else
				{
					zR.push_back(i < 0 ? - z[abs(i)] : z[i]);
					dPhidR_.push_back(dPhidR_makima[abs(i)](Rb));
					dPhidz_.push_back(dPhidz_makima[abs(i)](Rb));
				}
			}
			std::vector<double> zz = zR;
			boost::math::interpolators::makima dPhidR_spline(std::move(zR), std::move(dPhidR_));
			boost::math::interpolators::makima dPhidz_spline(std::move(zz), std::move(dPhidz_));
			f.fR = - dPhidR_spline(zb) * r_sq_ratio;
			f.fz = pos.z < 0 ? dPhidz_spline(zb) * r_sq_ratio: - dPhidz_spline(zb) * r_sq_ratio;
		}
		return f;
	}

	double CartGrid::pot(double r)
	{
		if(r < R_max) return Phi_makima[0](r);
		else return Phi_makima[0](R_max) * R_max / r;
	}

	double CartGrid::dPhidr(double r)
	{
		if(r < R_max) return dPhidR_makima[0](r);
		else return dPhidR_makima[0](R_max) * RC::sq(R_max / r);
	}

	double CartGrid::d2Phidr2(double r)
	{
		if(r < R_max) return dPhidR_makima[0].prime(r);
		else return dPhidR_makima[0].prime(R_max) * RC::cb(R_max / r);
	}
}
*/
