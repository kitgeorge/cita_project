/**
 * \file		shared/inc/potential_CartGrid.hpp
 * \brief		Cartesian grid.
 * \author	Rimpei Chiba
 * \date		2022-
 */
/*
#ifndef _RC_POTENTIAL_CARTGRID_H_
#define _RC_POTENTIAL_CARTGRID_H_

#include <stdlib.h>
#include <iostream>
#include <math.h>
#include <boost/math/interpolators/makima.hpp>

#include "math_RC.hpp"
#include "import.hpp"
#include "galaxy.hpp"
#include "force.hpp"
#include "coordinate.hpp"
#include "potential.hpp"
#include "NewtonRaphson.hpp"

namespace RC
{
  class CartGrid : public Potential
  {
    private:
      double R_max, z_max;
      std::vector<double> R{}, z{};
      std::vector<std::vector<double>> Phi{}, dPhidR{}, dPhidz{};
      std::vector<boost::math::interpolators::makima<std::vector<double>>> 
        Phi_makima{}, dPhidR_makima{}, dPhidz_makima{};
      void readAGAMA(const std::string &in, 
        std::vector<std::vector<double>> &Phi, 
        std::vector<std::vector<double>> &dPhidR, 
        std::vector<std::vector<double>> &dPhidz);

    public:
      double pot(const Polar &pos);
      Force force(const Polar &pos);
      double pot(double r);
      double dPhidr(double r);
      double d2Phidr2(double r);
      CartGrid(const std::string in);
  };
}

#endif
*/