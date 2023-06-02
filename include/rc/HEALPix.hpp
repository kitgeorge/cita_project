/**
 * \file		shared/src/HEALPix.hpp
 * \brief		Hierarchical Equal Area and isoLatitude Pixelization (HEALPix)
 *          Equations given in Gorski et al. (2005)
 * \author	Rimpei Chiba
 * \date		2020
 */

#ifndef _HEALPIX_H_
#define _HEALPIX_H_

#include <stdlib.h>
#include <vector>
#include <iostream>

#include "pi.hpp"
#include "utils.hpp"

struct sph_coord
{
  double r, theta, phi;
  sph_coord(double r, double theta, double phi):
  r(r), theta(theta), phi(phi)
  {};
};

class HEALPix
{
  public:
    const int N_base_theta, N_base_phi, N_base_pix,
              N_side, N_side_sq, N_pix, N_theta, N_phi_belt;
    const double base_pix_area, pix_area, dz, dz_cap, z_transit, phi_base, dzdphi;
    std::vector<int> N_phi, N_cap_pix;
    int N_cap_pix_tot;

    sph_coord get_pixel_centre(int i_pix);
    unsigned int get_pixel_index(double theta, double phi);
    unsigned int get_pixel_index_i(int i_theta, int i_phi);

    HEALPix(int N_base_theta, int N_base_phi, int N_side);
};

#endif
