/**
 * \file		shared/inc/NewtonRaphson.hpp
 * \brief		Find roots of equations using Newton-Raphson method.
 *          Extracted from the falcon package.
 * \author	Walter Dehnen
 * \date		2020
 */

#ifndef _NEWTONRAPSHON_H_
#define _NEWTONRAPSHON_H_

#include <stdlib.h>
#include <vector>
#include <math.h>
#include <ctype.h>

#include "utils.hpp"

template<typename FuncAndDeriv>
double BracketNewtonRaphson(double xl, double xh, FuncAndDeriv const& func, double xacc)
{
  const int maxit = 100;
  double df;
  auto fl = func(xl, df);
  auto fh = func(xh, df);
  if(fl * fh > 0)
  {
    show("Warning [shared/inc/NewtonRaphson.hpp]: root must be bracketed.");
    show("xl",xl);
    show("xh",xh);
    show("fl",fl);
    show("fh",fh);
  }
  if(fl > 0)
  {
    std::swap(xl, xh);
    std::swap(fl, fh);
  }
  auto rts = 0.5 * (xl + xh);
  auto dxo = std::abs(xh - xl);
  auto dx  = dxo;
  auto f   = func(rts, df);
  for(int j = 0; j != maxit; ++j)
  {
    if((((rts - xh) * df - f) * ((rts - xl) * df - f) >= 0.) ||
      (std::abs(f+f) > std::abs(dxo*df)))
    {
      dxo = dx;
      dx  = 0.5 * (xh - xl);
      rts = xl + dx;
      if(xl <= rts && xl >= rts) return rts;
    }
    else
    {
      // Newton-Raphson method
      dxo = dx;
      dx  = f / df;
      auto tmp = rts;
      rts -= dx;
      if(tmp <= rts && tmp >= rts) return rts;
    }
    // End iteration if dx is smaller than the tolerance
    if(std::abs(dx) < xacc) return rts;

    // Update either ends with the new root
    // (use this f and df for the next iteration)
    f = func(rts, df);
    if(f < 0)
    {
      xl = rts;
      fl = f;
    }
    else
    {
      xh = rts;
      fh = f;
    }
  }
  //show("BracketNewtonRaphson(): maximum number of " 
  // + std::to_string(maxit) +  " iterations exceeded.");
  return rts;
}

template<typename FuncAndDeriv>
double bisection(double xl, double xh, FuncAndDeriv const& func, double xacc)
{
  const int maxit = 100;
  auto fl = func(xl);
  auto fh = func(xh);
  if(fl * fh > 0)
  {
    show("Warning [shared/inc/NewtonRaphson.hpp/bisection]: root must be bracketed.");
    show("xl",xl);
    show("xh",xh);
    show("fl",fl);
    show("fh",fh);
  }
  if(fl > 0)
  {
    std::swap(xl, xh);
    std::swap(fl, fh);
  }
  double rts;
  for(int j = 0; j != maxit; ++j)
  {
    rts = 0.5 * (xl + xh);
    if(((double)func(rts) == 0.) || (std::abs(xh - xl) < xacc)) break;
    else if(func(rts) * func(xl) < 0) xh = rts;
    else xl = rts;
  }
  return rts;
}

#endif