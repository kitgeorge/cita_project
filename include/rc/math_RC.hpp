/**
 * \file		shared/inc/math_RC.hpp
 * \brief		Mathematic functions.
 * \author	Rimpei Chiba
 * \date		2020-
 */

#ifndef _MATHRC_H_
#define _MATHRC_H_

#include <stdlib.h>
#include <vector>
#include <limits>
#include <math.h>
#include <assert.h>
#include <complex>
#include <numeric>

#include "utils.hpp"

namespace RC
{
  template <typename T>
  double norm(std::vector<T> vec);
  template <typename T>
  T sq(const T x);
  template <typename T>
  T cb(const T x);
  template <typename T>
  T sgn(const T x); 

  unsigned factorial(unsigned n);
  unsigned doublefactorial(unsigned n);
  
  template<typename Real>
  Real gegenbauer(unsigned int n, Real a, Real x);
  template<typename Real>
  Real gegenbauer_derivative(unsigned int n, Real a, Real x, unsigned int k);

  double WignerdMatrix(unsigned l, int mp, int m, double beta);
  double WignerdMatrix(unsigned l, int mp, int m, double beta, double &derivative);
  std::complex<double> WignerDMatrix(unsigned l, int mp, int m, 
    double alpha, double beta, double gamma);
  double WignerdMatrix22(int mp, double beta);
  double WignerdMatrix22(int mp, double beta, double &dddcosb);
  
  double discreteDelta(double w, double x);
  double HeavisideStep(double x);
}

#endif