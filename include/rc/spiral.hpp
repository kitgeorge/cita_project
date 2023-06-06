/**
 * \file		shared/inc/spiral.hpp
 * \brief		Functions for spiral model.
 * \author	Rimpei Chiba
 * \date		2021
 */
/*
#ifndef _SPIRAL_H_
#define _SPIRAL_H_

#include <stdlib.h>
#include <math.h>
#include <boost/math/special_functions/spherical_harmonic.hpp>

#include "galaxy.hpp"
#include "force.hpp"
#include "utils.hpp"

// Model of Spiral
// 1. tightly wound approximation
// 2. rigid rotation
// 3. transient

namespace RC
{
  class Spiral
  {
    public:
      int m;
      double A, Rs, pitch, phi0, Omegap, t0, ts_sq_2;
    
    public:
      virtual double pot(double R, double phi, double t);
      virtual RC::Force force(double R, double phi, double t);

      Spiral(int m, double A, double Rs, double pitch, 
        double phi0, double Omegap, double t0, double ts);
  };

  class Spiral2 : Spiral
  {
    public:
      double pot(double R, double phi, double t);
      RC::Force force(double R, double phi, double t);

      Spiral2(int m, double A, double Rs, double pitch, 
        double phi0, double Omegap, double t0, double ts) : 
        Spiral(m, A, Rs, pitch, phi0, Omegap, t0, ts)
      {
      }
  };
}

#endif
*/