/**
 * \file		shared/inc/force.cpp
 * \brief		Force.
 * \author	Rimpei Chiba
 * \date		2020-
 */

#include "force.hpp"

namespace RC
{
  void Force::amp(double x)
  {
    fR   *= x;
    fphi *= x;
    fz   *= x;
  }
  void Force::add(Force F)
  {
    fR   += F.fR;
    fphi += F.fphi;
    fz   += F.fz;
  }
}

