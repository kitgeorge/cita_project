/**
 * \file		shared/src/sht.cpp
 * \brief		Spherical harmonics transform. Not optimized. 
 * \author	Rimpei Chiba
 * \date		2020 
 */

#include "sht.hpp"

//
// map f(θ,φ) to alm
//
// f(θ,φ) = Σ_lm alm Ylm(θ,φ)
// alm = ∫∫ dθdφsinθ f(θ,φ) Ylm(θ,φ)^*  (* conjugate)
//     = ∫∫ dθdφsinθ f(θ,φ) Plm(θ) exp(-imφ)
//     = ∫ dθsinθ Plm(θ) ∫ dφ f(θ,φ) exp(-imφ)
//     = ∫ dθsinθ Plm(θ) Fm(θ)   (FFT)
//
// Example code given at the end
//

void map2alm_ECP(unsigned int l_max, unsigned int m_max, 
  std::vector<std::vector<std::complex<double>>> &f, 
  std::vector<std::vector<std::complex<double>>> &alm)
{
  alm.resize(l_max);
  std::fill(alm.begin(), alm.end(), 
    std::vector<std::complex<double>>(m_max, std::complex<double>(0.,0.)));
  const auto N_theta = f.size();
  const auto dtheta = Pi / N_theta;
  for(size_t j = 0; j < N_theta; ++j)
  {
    const auto theta = dtheta * (j + 0.5);
    const auto cos_theta = cos(theta);
    const auto sin_theta = sin(theta);

    // FFT in azimuth phi
    const auto N_phi = f[j].size();
    const auto dphi = TPi / N_phi;
    fftw_complex *in, *out;
    fftw_plan p;
    in  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N_phi);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N_phi);
    for(size_t k = 0; k < N_phi; ++k)
    {
      in[k][0] = f[j][k].real();
      in[k][1] = f[j][k].imag();
    }
    p = fftw_plan_dft_1d(N_phi, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(p);

    const size_t M = m_max < N_phi ? m_max : N_phi;
    for(size_t m = 0; m < M; ++m)
    {
      std::complex<double> Fm(out[m][0], out[m][1]);
      for(size_t l = m; l < l_max; ++l)
      {
        double factorial_neg = 1.;
        double factorial_pos = 1.;
        for(size_t i = 1; i <= l-m; ++i) factorial_neg *= i;
        for(size_t i = 1; i <= l+m; ++i) factorial_pos *= i;
        const auto coeff = sqrt((l+l+1)*iFPi*factorial_neg/factorial_pos);
        const auto Plm = boost::math::legendre_p(l, m, cos_theta);
        alm[l][m] += coeff * Plm * Fm * sin_theta * dtheta * dphi;
      }
    }
  }
}

void map2alm(unsigned int l_max, unsigned int m_max, 
  std::vector<std::vector<std::complex<double>>> &f, std::vector<double> &theta_v, 
  std::vector<std::vector<std::complex<double>>> &alm)
{
  alm.resize(l_max);
  std::fill(alm.begin(), alm.end(), 
    std::vector<std::complex<double>>(m_max, std::complex<double>(0.,0.)));
  
  const auto N_theta = f.size();
  for(size_t j = 0; j < N_theta; ++j)
  {
    // FFT in azimuth phi
    const auto N_phi = f[j].size();
    const auto dphi = TPi / N_phi;
    fftw_complex *in, *out;
    fftw_plan p;
    in  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N_phi);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N_phi);
    for(size_t k = 0; k < N_phi; ++k)
    {
      in[k][0] = f[j][k].real();
      in[k][1] = f[j][k].imag();
    }
    p = fftw_plan_dft_1d(N_phi, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(p);

    // Get theta,dtheta from theta index j
    const auto theta = theta_v[j];
    const auto cos_theta = cos(theta);
    const auto sin_theta = sin(theta);
    const auto dtheta = j == 0 ? 0.5 * (theta_v[1] + theta_v[0]) : 
                        j == N_theta - 1 ? Pi - 0.5 * (theta_v[j] + theta_v[j-1]) :
                        0.5 * (theta_v[j+1] - theta_v[j-1]);
    const int M = m_max < N_phi ? m_max : N_phi;
    for(int m = 0; m < M; ++m)
    {
      std::complex<double> Fm(out[m][0], out[m][1]);
      for(size_t l = m; l < l_max; ++l)
      {
        double factorial_neg = 1.;
        double factorial_pos = 1.;
        for(size_t i = 1; i <= l-m; ++i) factorial_neg *= i;
        for(size_t i = 1; i <= l+m; ++i) factorial_pos *= i;
        const auto coeff = sqrt((l+l+1)*iFPi*factorial_neg/factorial_pos);
        const auto Plm = boost::math::legendre_p(l, m, cos_theta);
        alm[l][m] += coeff * Plm * Fm * sin_theta * dtheta * dphi;
      }
    }
  }
}

//
// Test code
//
/*
int main(int argc, char ** argv)
{
  if(argc != 3) std::cout << "Provide l and m: ex.) ./sht 2 2" << std::endl;
  else
  {
    // Make f (any function on a sphere)
    // Ex.) spherical harmonics of mode l,m specified by the command line arguments
    std::vector<std::vector<std::complex<double>>> f;
    std::vector<double> theta_v;
    const int l = atoi(argv[1]);
    const int m = atoi(argv[2]);
    const int N_theta = 100;
    const int N_phi   = 100;
    const auto dtheta = Pi / N_theta;
    const auto dphi   = TPi / N_phi;

    double factorial_neg = 1.;
    double factorial_pos = 1.;
    for(int i = 1; i <= l-m; ++i) factorial_neg *= i;
    for(int i = 1; i <= l+m; ++i) factorial_pos *= i;
    const auto coeff = sqrt((l+l+1)*iFPi*factorial_neg/factorial_pos);

    for(int i = 0; i < N_theta; ++i)
    {
      const auto theta = dtheta * (i + 0.5);
      theta_v.push_back(theta);
      const auto Plm = boost::math::legendre_p(l, m, cos(theta));
      std::vector<std::complex<double>> buf;
      for(int j = 0; j < N_phi; ++j)
      {
        const auto phi = dphi * j;
        std::complex<double> exp_i_m_phi(cos(m * phi), sin(m * phi));
        buf.push_back(coeff * Plm * exp_i_m_phi);
      }
      f.push_back(buf);
    }

    // Map to a_lm
    unsigned int l_max = 7;
    unsigned int m_max = 7; 
    std::vector<std::vector<std::complex<double>>> alm;
    map2alm(l_max, m_max, f, theta_v, alm);

    // Export results
    for(size_t l = 0; l < alm.size(); ++l)
    {
      for(size_t m = 0; m < alm[l].size(); ++m)
      {
        if(fabs(alm[l][m].real()) < 1.e-8) alm[l][m] = std::complex<double>(0.,alm[l][m].imag());
        if(fabs(alm[l][m].imag()) < 1.e-8) alm[l][m] = std::complex<double>(alm[l][m].real(),0.);
        std::cout << "a[" << l << "][" << m << "] = " << alm[l][m] << std::endl;
      }
      std::cout << std::endl;
    }
  }
  return EXIT_SUCCESS;
}
*/
