/**
 * \file		shared/inc/force.hpp
 * \brief		Force.
 * \author	Rimpei Chiba
 * \date		2020-
 */

#ifndef _FORCE_H_
#define _FORCE_H_

#include <stdlib.h>

namespace RC
{
  struct Force
  {
    double fR{}, fphi{}, fz{};
    void amp(double x);
    void add(Force F);

    Force operator + (const Force &F)
    {
      Force temp;
      temp.fR   = fR   + F.fR;
      temp.fphi = fphi + F.fphi;
      temp.fz   = fz   + F.fz;
      return temp;
    }
    void operator += (const Force &F)
    {
      fR   += F.fR;
      fphi += F.fphi;
      fz   += F.fz;
    }
    void operator -= (const Force &F)
    {
      fR   -= F.fR;
      fphi -= F.fphi;
      fz   -= F.fz;
    }
    void operator *= (const double &amp)
    {
      fR   *= amp;
      fphi *= amp;
      fz   *= amp;
    }
    void operator /= (const double &amp)
    {
      fR   /= amp;
      fphi /= amp;
      fz   /= amp;
    }
  };
}

#endif