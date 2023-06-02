/**
 * \file galaxy.hpp
 * \brief Properties of the Sun/Galaxy.
 * 
 */

#ifndef _GALAXY_H_
#define _GALAXY_H_

#include "units.hpp"

const double R_sun     = 8.2 * Units::kpc;            // Galactocentric radius of the Sun. McMillan 2017.
const double R_sun_sq  = R_sun * R_sun;
const double phi_sun   = - 30. * Units::degree;       // Azimuthal angle of the Sun w.r.t the bar. Wegg 2015.
const double z_sun     = 0.02 * Units::kpc;           // Vertical height of the Sun. Joshi 2007.
const double vR_sun    = -11.1 * Units::kms;          // Radial velocity of the Sun. Schoenrich 2010.
const double vc_sun    = 235. * Units::kms;           // Circular velocity at R_sun. Reid 2019.
const double vc_sun_sq = vc_sun * vc_sun;
const double vphi_sun  = 12.24 * Units::kms + vc_sun; // Azimuthal velocity of the Sun. Schoenrich 2010.
const double vz_sun    = 7.25 * Units::kms;           // Vertical velocity of the Sun. Schoenrich 2010.
const double Omega_sun = vc_sun / R_sun;              // Circular frequency at R_sun.

#endif
