/**
 * \file		shared/src/HEALPix.cpp
 * \brief		Hierarchical Equal Area and isoLatitude Pixelization (HEALPix)
 *          Equations given in Gorski et al. (2005)
 * \author	Rimpei Chiba
 * \date		2020 
 */

#include "HEALPix.hpp"

HEALPix::HEALPix(int N_base_theta, int N_base_phi, int N_side):
  N_base_theta(N_base_theta),          // Number of base pixels in theta
  N_base_phi(N_base_phi),              // Number of base pixels in phi
  N_base_pix(N_base_theta*N_base_phi), // Total number of base pixels
  N_side(N_side),                      // Number of division of the base pixels on one side
  N_side_sq(N_side * N_side),          // Number of pixels in a single base pixel
  N_pix(N_base_pix * N_side_sq),       // Total number of pixels
  N_theta((N_base_theta+1)*N_side-1),  // Number of isolatitudinal rings
  N_phi_belt(N_base_phi * N_side),     // Number of azimuthal pixels at the belt
  base_pix_area(FPi / N_base_pix),     // Area of base pixel
  pix_area(base_pix_area / N_side_sq), // Area of pixel
  dz(4. / (N_base_theta * N_side)),    // Vertical width of pixel (rhombus)
  dz_cap(1. / N_base_theta),           // Vertical width of the 'cap'
  z_transit(1. - dz_cap),              // Transit z between polar 'cap' and equatorial 'belt'
  phi_base(TPi / N_base_phi),          // Azimuthal width of base pixels
  dzdphi(4./(N_base_theta*phi_base)),  // Inclination of the rhombus
  // Vectors. Redefine in constructor.
  N_phi(N_theta, N_phi_belt),          // Number of azimuthal pixels at each latitude
  N_cap_pix(N_side - 1, 0)             // Number of pixels in the polar 'cap'
{
  for(int k = 0; k < N_side - 1; ++k)
  {
    N_phi[k] = N_base_phi * (k + 1);
    N_phi[N_theta - 1 - k] = N_phi[k];
    if(k == 0) N_cap_pix[k] = N_phi[k];
    else if (k > 0) N_cap_pix[k] = N_cap_pix[k-1] + N_phi[k];
  }
  N_cap_pix_tot = N_cap_pix[N_side - 2];
}

sph_coord HEALPix::get_pixel_centre(int i_pix)
{
  if((i_pix < 0) || (N_pix <= i_pix))
  {
    std::cout << "Warning [shared/src/HEALPix.cpp]: Pixel index i_pix=" 
      << i_pix << " is out of range [0," << N_pix << ")" << std::endl;
    sph_coord sph(0,0,0);
    return sph;
  }
  else
  {
    // i_pix = pixel index in RING scheme (0 <= i_pix < N_pix)
    int i{}; // latitudinal index (0 <= i < N_theta)
    int j{}; // azimuthal index (0 <= j < N_phi[i])
    double z{}, phi{};

    if(i_pix < N_cap_pix_tot) // pixel in north cap
    {
      // Find i,j from i_pix
      for(int k = 0; k < N_side - 1; ++k)
      {
        if(i_pix < N_cap_pix[k])
        {
          i = k;
          if(k == 0) j = i_pix;
          else j = i_pix - N_cap_pix[k-1];
          break;
        }
      }
      // Get theta,phi from i,j
      z = 1. - (double)(i+1.)*(i+1.)*dz_cap/N_side_sq;
      phi = TPi / N_phi[i] * (j + 0.5);
    }
    else if(i_pix < N_pix - N_cap_pix_tot) // pixel in north/south belt
    {
      // Find i,j from i_pix
      i = N_side - 1 + floor((i_pix - N_cap_pix_tot) / N_phi_belt);
      j = (i_pix - N_cap_pix_tot) % N_phi_belt;
      // Get theta,phi from i,j
      z = 1. + dz_cap * (1. - 2.*(i+1.)/N_side);
      auto s = (i - N_side + 2) % 2; // (i + 1 - (N_side - 1)) mod 2
      phi = TPi / N_phi[i] * (j + 0.5 * s);
    }
    else if(i_pix < N_pix) // pixel in south cap
    {
      int i_pix_back = N_pix - i_pix - 1;
      // Find i,j from i_pix_back
      for(int k = 0; k < N_side - 1; ++k)
      {
        if(i_pix_back < N_cap_pix[k])
        {
          i = N_theta - k - 1;
          if(k == 0) j = N_phi[k] - i_pix_back - 1;
          else j = N_phi[k] - (i_pix_back - N_cap_pix[k-1]) - 1;
          break;
        }
      }
      // Get theta,phi from i,j
      z = - 1. + (double)(N_theta-i)*(N_theta-i)*dz_cap/N_side_sq;
      phi = TPi / N_phi[i] * (j + 0.5);
    }
    const auto theta = acos(z);
    sph_coord sph(0, theta, phi);
    return sph;
  }
}

unsigned int HEALPix::get_pixel_index(double theta, double phi)
{
  const auto z = cos(theta);
  phi = Oto2PI(phi);
  const int i_phi_base = floor(phi / phi_base); // quotient
  const auto phit = std::fmod(phi, phi_base); // remainder
  if(z >= z_transit) // North cap
  {
    int k, l;
    // upper-left boundary
    for(k = 0; k < N_side - 1; ++k)
    {
      auto buf = (k + 1.) * phi_base / (N_side * phit);
      auto zb = 1. - buf * buf * dz_cap;
      if(z >= zb) break;
    }
    // upper-right boundary
    for(l = 0; l < N_side - 1; ++l)
    {
      auto buf = (l + 1.) * phi_base / (N_side * (phit - phi_base));
      auto zb = 1. - buf * buf * dz_cap;
      if(z >= zb) break;
    }
    // Find i_pix from k,l
    const auto m = l + k;
    int i_pix{};
    for(int i = 0; i <= m; ++i) i_pix += i;
    i_pix = (N_base_phi * i_pix) + k + i_phi_base * (m + 1);
    return i_pix;
  }
  else if (z >= - z_transit) // North/south belt
  {
    // upper-left boundary
    int k = 0;
    while(z < z_transit - dz*(k+1.) + dzdphi*phit) ++k;
    // upper-right boundary
    int l = 0;
    while(z < z_transit - dz*(l+1.) - dzdphi*(phit-phi_base)) ++l;
    // Find i_pix from k,l
    const auto n = l + k - N_side + 1;
    int i_pix = N_cap_pix_tot + i_phi_base * N_side 
      + n * N_side * N_base_phi - floor(n * 0.5) + k;
    if((i_phi_base == N_base_phi - 1) && (k - l == N_side))
      i_pix -= N_phi_belt; // edge at phi=2pi
    return i_pix;
  }
  else // South cap
  {
    int k, l;
    // upper-left boundary
    for(k = 0; k < N_side - 1; ++k)
    {
      auto buf = (k + 1.) * phi_base / (N_side * (phit - phi_base));
      auto zb = - 1. + buf * buf * dz_cap;
      if(z <= zb) break;
    }
    // upper-right boundary
    for(l = 0; l < N_side - 1; ++l)
    {
      auto buf = (l + 1.) * phi_base / (N_side * phit);
      auto zb = - 1. + buf * buf * dz_cap;
      if(z <= zb) break;
    }
    // Find i_pix from k,l
    const auto m = l + k;
    int i_pix{};
    for(int i = 0; i <= m; ++i) i_pix += i;
    i_pix = N_pix - 1 - ((N_base_phi * i_pix) + k 
      + (N_base_phi - 1 - i_phi_base) * (m + 1));
    return i_pix;
  }
}

unsigned int HEALPix::get_pixel_index_i(int i_theta, int i_phi)
{
  unsigned int i_pix = 0;
  if(i_theta >= N_theta)
  {
    std::cout << "Warning [shared/src/HEALPix.cpp]: theta index i_theta=" 
      << i_theta << " is out of range [0," << N_theta << ")" << std::endl;
  }
  else if(i_phi >= N_phi[i_theta])
  {
    std::cout << "Warning [shared/src/HEALPix.cpp]: phi index i_phi=" 
      << i_phi << " is out of range [0," << N_phi[i_theta] << ")" << std::endl;
  }
  else
  {
    if(i_theta < N_side - 1)
    {
      for(int k = 0; k < i_theta; ++k) i_pix += k + 1;
      i_pix *= 4;
      i_pix += i_phi;
    }
    else if(i_theta < N_theta - (N_side - 1))
    {
      i_pix = N_cap_pix_tot + (i_theta - (N_side - 1)) * N_phi_belt + i_phi;
    }
    else
    {
      const auto i_theta_rev = N_theta - i_theta - 1;
      int i_pix_reverse = 0;
      for(int k = 0; k < i_theta_rev; ++k) i_pix_reverse += k + 1;
      i_pix_reverse *= 4;
      i_pix_reverse += N_phi[i_theta] - i_phi;
      i_pix = N_pix - i_pix_reverse;
    }
  }
  return i_pix;
}

/*
//
// Test code
//
int main()
{
  const int N_base_theta = 3;
  const int N_base_phi = 4;
  const int N_side = 4;
  HEALPix healpix(N_base_theta, N_base_phi, N_side);

  for(int i = 0; i < healpix.N_pix; ++i)
  {
    sph_coord sph = healpix.get_pixel_centre(i);
    unsigned int i_out = healpix.get_pixel_index(sph.theta, sph.phi);
    double z = cos(sph.theta);
    // Confirm i_in == i_out
    std::cout << "\ni_in   = " << i 
              << "\nphi/Pi = " << sph.phi / Pi
              << "\nz      = " << z
              << "\ni_out  = " << i_out << std::endl;
  }

  for(int i = 0; i < healpix.N_theta; ++i)
  {
    for(int j = 0; j < healpix.N_phi[i]; ++j)
    {
      unsigned int i_pix = healpix.get_pixel_index_i(i,j);
      show("i_pix",i_pix);
    }
    std::cout << std::endl;
  }
  return EXIT_SUCCESS;
}
*/
