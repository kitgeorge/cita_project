/**
 * \file		shared/src/sht.hpp
 * \brief		Spherical harmonics transform. Not optimized.
 * \author	Rimpei Chiba
 * \date		2020
 */

#ifndef _SPHERICAL_HARMONICS_TRANSFORM_H_
#define _SPHERICAL_HARMONICS_TRANSFORM_H_

#include <stdlib.h>
#include <vector>
#include <cmath>
#include <complex>
#include <boost/math/special_functions/spherical_harmonic.hpp>
#include <boost/math/special_functions/legendre.hpp>
#include <fftw3.h>

#include "pi.hpp"

// Works with any isolatitude sampling methods
// e.g. Equidistant Cylindrical Projection (ECP),
// Hierarchical Equal Area and isoLatitude Pixelization (HEALPix)

void map2alm(unsigned int l_max, unsigned int m_max, 
  std::vector<std::vector<std::complex<double>>> &f, std::vector<double> &theta_v, 
  std::vector<std::vector<std::complex<double>>> &alm);

#endif