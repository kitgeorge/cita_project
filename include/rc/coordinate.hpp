/**
 * \file		shared/inc/coordinate.hpp
 * \brief		Coordinates
 * \author	Rimpei Chiba
 * \date		2022-
 */

#ifndef _COORDINATE_H_
#define _COORDINATE_H_

#include <stdlib.h>

namespace RC
{
  class Polar
  {
    public:
      double x{}, y{}, z{}, R{}, phi{}, r{}, theta{},
        cp{}, sp{}, ct{}, st{}, cot{};
      void rot(double p) { phi -= p; };

      Polar(double x, double y, double z) : x(x), y(y), z(z)
      {
        R = x * x + y * y;
        r = sqrt(R + z * z);
        R = sqrt(R);
        phi = atan2(y, x);
        theta = asin(R / r);
        cp = x / R; // cos(φ)
        sp = y / R; // sin(φ)
        ct = z / r; // cos(θ)
        st = R / r; // sin(θ)
        cot = z / R; // cot(θ)
      };
  };
}

#endif