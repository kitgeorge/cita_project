/**
 * \file		shared/inc/bar.cpp
 * \brief		Functions for bar models.
 * \author	Rimpei Chiba
 * \date		2020-
 */

#include "bar.hpp"

namespace RC
{
  double BarGrowthDehnen(double t, double Tg)
  {
    if(t < 0) return 0.;
    else if(t < Tg)
    {
      const auto xi = (t + t) / Tg - 1.;
      const auto xi_2 = xi * xi;
      const auto xi_3 = xi_2 * xi;
      const auto xi_5 = xi_3 * xi_2;
      // 3 / 16  = 0.1875
      // 5 / 8   = 0.625
      // 15 / 16 = 0.9375
      return 0.1875 * xi_5 - 0.625 * xi_3 + 0.9375 * xi + 0.5;
    }
    else return 1.;
  }

  double BarGrowthCosine(double t, double Tg)
  {
    if(t < 0) return 0.;
    else if(t < Tg) return 0.5 * (1. - cos(Pi * t / Tg));
    else return 1.;
  }

  // ----- Bar base class -----
  double Bar::getwb(double t)
  {
    if(eta == 0) return wb0;
    else
    {
      if(t < Tg) return wb0;
      else if(t < t2)
      {
        double xi = (t - Tg) / Tt;
        return wb0 / (1. + C_tr_sq * xi * xi);
      }
      // General solution to eta = - dwb/dt 1/wb^2 = const.
      // wb = wb0 / (1 + eta wb0 t).
      else return wb0 / (1. + C_tr_sq + wb0_eta * (t - t2));
    }
  }

  double Bar::getdwbdt(double t)
  {
    if(eta == 0) return 0.;
    else
    {
      if(t < Tg) return 0.;
      else if(t < t2) return - eta * RC::sq(getwb(t)) * (t - Tg) / Tt;
      else return - eta * RC::sq(getwb(t));
    }
  }

  double Bar::getphib(double t)
  {
    // Integrate wb(t)
    if(t < Tg) return phib0 + wb0 * t;
    else
    {
      if(eta == 0) return phib0 + wb0 * t;
      else if(Tt == 0) return phib0 + wb0 * Tg + log(1. + wb0_eta * (t - Tg)) / eta;
      else
      {
        if(t < t2) return phib0 + wb0 * (Tg + Tt * C_tr_i * atan(C_tr * (t - Tg) / Tt));
        else return phib0 + wb0 * (Tg + Tt * C_tr_i * atan(C_tr) 
          + log(1. + wb0_eta * (t - t2) / (1. + C_tr_sq)) / wb0_eta);
      }
    }
  }




  void Bar::setphib0(double phib, double t)
  {
    if(t < Tg) phib0 = phib - wb0 * t;
    else
    {
      if(eta == 0) phib0 = phib - wb0 * t;
      else if(Tt == 0) phib0 = phib - (wb0 * Tg + log(1. + wb0_eta * (t - Tg)) / eta);
      else
      {
        if(t < t2) phib0 = phib - wb0 * (Tg + Tt * C_tr_i * atan(C_tr * (t - Tg) / Tt));
        else phib0 = phib - wb0 * (Tg + Tt * C_tr_i * atan(C_tr) 
          + log(1. + wb0_eta * (t - t2) / (1. + C_tr_sq)) / wb0_eta);
      }
    }
  }

  double Bar::Phim(double r)
  {
    (void)r;
    return 0.;
  }

  double Bar::dPhimdr(double r)
  {
    (void)r;
    return 0.;
  }

  double Bar::Philm(double r)
  {
    (void)r;
    return 0.;
  }

  double Bar::density(const Polar &pos)
  {
    (void)pos;
    return 0.;
  }

  double Bar::pot(const Polar &pos)
  {
    (void)pos;
    return 0.;
  }

  RC::Force Bar::force(const Polar &pos)
  {
    (void)pos;
    Force f;
    return f;
  }

  void Bar::update(double new_wb, double new_rCR)
  {
    wb = new_wb;
    rCR = new_rCR;
  }

  std::vector<std::vector<std::vector<std::complex<double>>>>
    Bar::calcblmn(int rank, int size)
  {
    (void) rank;
    (void) size;
    return blmn;
  }

  void Bar::setEvolPara(double Tg_, double Tt_, double wb0_, double eta_)
  {
    Tg = Tg_;
    Tt = Tt_;
    wb0 = wb0_;
    eta = eta_;
    t2 = Tg + Tt;
    wb0_eta = wb0 * eta;
    C_tr_sq = 0.5 * Tt * wb0_eta;
    C_tr = sqrt(C_tr_sq);
    C_tr_i = 1. / C_tr;
  }


  //
  // ----- Y22 Bar -----
  // (Chiba et al. 2020, Chiba & Schoenrich 2021)
  //
  // Phib(x) = Phim(r) * (sinθ)^2 * cos(m * φRF)
  //
  // φRF     = φ - φb
  // Phim(r) = Phim_C * a^2 * ((b+1)/(b+a))^5
  // Phim_C  = - A vc_sun^2 / m
  // a       = r / rCR
  // b       = rb / rCR
  //
  // dPhim/dr = da/dr dPhim/da
  //          = 1/rCR * Phim_C * a * ((b+1)/(b+a))^5 * (2b-3a)/(b+a)
  //          = Phim_C / r * a^2 * ((b+1)/(b+a))^5 * (2b-3a)/(b+a)
  //          = Phim(r) / r * (2b-3a)/(b+a)
  //
  // fr = - dPhib/dr
  //    = - dPhim/dr * (sinθ)^2 * cos(m * φRF)
  // fθ = - 1/r dPhib/dθ
  //    = - 1/r * Phim(r) * 2*sinθcosθ * cos(m * φRF)
  // fφ = - 1/(r sinθ) dPhib/dφ
  //    = - 1/R dPhib/dφ
  //    = 1/R * m * Phim(r) * (sinθ)^2 * sin(m * φRF)
  //

  Y22Bar::Y22Bar(int m, double A, double b):
    m(m),
    b(b),
    Phim_C(- A * vc_sun_sq / m)
  {}

  double Y22Bar::Phim(double r)
  {
    const auto a = r / rCR;
    return Phim_C * a * a * pow((b + 1.) / (b + a), 5);
  }

  double Y22Bar::dPhimdr(double r)
  {
    const auto a = r / rCR;
    return Phim(r) / r * ((b+b)-(a+a+a)) / (b+a);
  }

  double Y22Bar::Philm(double r)
  {
    const auto c = 1.2944172750371328; // 1 / (2 * Ylm(Pi/2,0))
    return Phim(r) * c;
  }

  double Y22Bar::density(const Polar &pos)
  {
    const auto rb = rCR * b;
    const auto x = pos.r / rb;
    const auto rhob = - Phim_C * b * b * pow(1. + 1. / b, 5)
      / (4. * Pi * Units::G * rb * rb) * 30. * x / pow(1. + x, 7);
    return rhob * RC::sq(pos.st) * cos(m * pos.phi);
  }

  double Y22Bar::pot(const Polar &pos)
  { 
    return Phim(pos.r) * RC::sq(pos.st) * cos(m * pos.phi);
  }

  RC::Force Y22Bar::force(const Polar &pos)
  {
    const auto a = pos.r / rCR;
    const auto _Phim = Phim(pos.r);
    const auto dPhimdR = _Phim / pos.r * ((b+b)-(a+a+a)) / (b+a);
    const auto m_phi = m * pos.phi;
    const auto st_sq = pos.st * pos.st;
    const auto fr = - dPhimdR * st_sq * cos(m_phi);
    const auto ftheta = (pos.z == 0) ? 0. : - _Phim * (pos.st + pos.st) * pos.ct * cos(m_phi) / pos.r;
    RC::Force f;
    f.fR   = fr * pos.st + ftheta * pos.ct;
    f.fz   = fr * pos.ct - ftheta * pos.st;
    f.fphi = m * _Phim * st_sq * sin(m_phi) / pos.R;
    return f;
  }

  //
  // ----- Ferrers Bar -----
  // (Chandrasekahr 1969, Ellipsoidal Figures of Equilibrium)
  // (Used in Chiba 2022)
  //
  FerrersBar::FerrersBar(int p, double a, double b, double c, double Mb, double rCR,
    int lmax, int nmax, double rs, double M):
    p(p),
    rCR_a(rCR / a),
    b_a(b / a),
    c_a(c / a),
    rhob0(Mb*(p+1)*doublefactorial(p+p+3)/(2.*doublefactorial(p+p+2)*Pi*a*b*c)),
    C_Ib(Mb/((p+p+5)*a*b*c)),
    a(a),
    b(b),
    c(c),
    lmax((lmax < 2) ? 2 : (lmax % 2 == 0) ? lmax : lmax - 1),
    nmax(nmax),
    N_grid(40),
    rs(rs),
    M(M),
    vgsq(Units::G * M / rs),
    Iln_i(lmax + 1, std::vector<double>(nmax + 1, 0.))
  {
    // Coefficients of basis functions
    blmn.resize(lmax + 1, std::vector<std::vector<std::complex<double>>>
    (lmax + lmax + 1, std::vector<std::complex<double>>(nmax + 1, std::complex<double>{})));
    // Iln for Hernquist & Ostriker radial basis function 
    for(int l = 2; l <= lmax; l += 2)
    {
      const auto l2 = l + l;
      for(int n = 0; n <= nmax; ++n)
      {
        const auto Kln = 0.5 * n * (n + l2 + l2 + 3) + (l + 1) * (l2 + 1);
        const auto g = tgamma(l2 + 1.5);
        const auto Iln = - Kln * FPi / pow(2, 8 * l + 6)
          * tgamma(n + l2 + l2 + 3) / (RC::factorial(n) * (n + l2 + 1.5) * g * g);
        Iln_i[l][n] = 1. / Iln;
      }
    }
  }

  double FerrersBar::densityxyz(double x, double y, double z)
  {
    const auto musq = sq(x/a) + sq(y/b) + sq(z/c);
    return musq < 1 ? rhob0 * pow(1. - musq, p) : 0;
  }

  double FerrersBar::density(const Polar &pos)
  {
    const auto musq = sq(pos.x/a) + sq(pos.y/b) + sq(pos.z/c);
    return musq < 1 ? rhob0 * pow(1. - musq, p) : 0;
  }

  double FerrersBar::pot(const Polar &pos)
  {
    const auto rn = pos.r / rs;
    const auto xi = (rn - 1.) / (rn + 1.);
    double pot{};
    for(int l = 2; l <= lmax; l += 2)
    {
      std::vector<double> Uln(nmax + 1, 0.);
      for(int n = 0; n <= nmax; ++n)
      {
        Uln[n] = - TSPi * pow(rn, l) / pow(1. + rn, l + l + 1)
          * RC::gegenbauer(n, l + l + 1.5, xi);
      }
      for(int m = - l; m <= l; m += 2)
      {
        if(m != 0)
        {
          const auto Ylm = boost::math::spherical_harmonic(l, m, pos.theta, pos.phi);
          for(int n = 0; n <= nmax; ++n)
          {
            double buf = Uln[n] * (blmn[l][m+lmax][n] * Ylm).real();
            if(!isnan(buf) && !isinf(buf)) pot += buf;
          }
        }
      }
    }
    return pot * Units::G * M / rs;
  }

  RC::Force FerrersBar::force(const Polar &pos)
  {
    const auto rn = pos.r / rs;
    const auto xi = (rn - 1.) / (rn + 1.);
    RC::Force f{};
    double fr{}, ft{}; // f_r, f_θ
    for(int l = 2; l <= lmax; l += 2)
    {
      std::vector<double> Uln(nmax + 1, 0.);
      std::vector<double> dUlndr(nmax + 1, 0.);
      const auto fl = (l == 0) ? - 1. : pow(rn, l - 1) * (1. - 2. * l * rn);
      for(int n = 0; n <= nmax; ++n)
      {
        const auto C = RC::gegenbauer(n, l + l + 1.5, xi);
        Uln[n]    = - TSPi * pow(rn, l) / pow(1. + rn, l + l + 1) * C;
        dUlndr[n] = - TSPi / (rs * pow(1. + rn, 2 * (l + 1))) * (fl * C
          + 2 * pow(rn, l) / (1. + rn) * RC::gegenbauer_derivative(n, l + l + 1.5, xi, 1));
      }
      for(int m = - l; m <= l; m += 2)
      {
        if(m != 0)
        {
          const auto Ylm  = boost::math::spherical_harmonic(l, m, pos.theta, pos.phi);
          const auto Ylm1 = boost::math::spherical_harmonic(l, m + 1, pos.theta, pos.phi);
          const auto dYlmdtheta = m * pos.cot * Ylm + sqrt((l - m) * (l + m + 1))
            * std::exp(std::complex<double>(0., - pos.phi)) * Ylm1;
          const auto dYlmdphi = std::complex<double>(0., m) * Ylm;
          for(int n = 0; n <= nmax; ++n)
          {
            const auto fr_ = - dUlndr[n] * (blmn[l][m+lmax][n] * Ylm).real();
            const auto ft_ = - Uln[n] * (blmn[l][m+lmax][n] * dYlmdtheta).real();
            const auto fp_ = - Uln[n] * (blmn[l][m+lmax][n] * dYlmdphi).real();
            if(!isnan(fr_) && !isinf(fr_)) fr     += fr_;
            if(!isnan(ft_) && !isinf(ft_)) ft     += ft_;
            if(!isnan(fp_) && !isinf(fp_)) f.fphi += fp_;
          }
        }
      }
    }
    ft /= pos.r;
    f.fphi /= pos.R;
    if(pos.z == 0) ft = 0;
    f.fR = fr * pos.st + ft * pos.ct;
    f.fz = fr * pos.ct - ft * pos.st;
    f.amp(vgsq);
    return f;
  }

  void FerrersBar::update(double new_wb, double new_rCR)
  {
    wb = new_wb;
    rCR = new_rCR;
    a = rCR / rCR_a;
    b = a * b_a;
    c = a * c_a;
    const auto C = C_Ib * a * b * c;
    Mb = C * (p + p + 5);
    Ib = C * (sq(a) + sq(b));
  }

  std::vector<std::vector<std::vector<std::complex<double>>>>
    FerrersBar::calcblmn(int rank, int size)
  {
    std::fill(blmn.begin(), blmn.end(), std::vector<std::vector<std::complex<double>>>
      (blmn[0].size(), std::vector<std::complex<double>>
      (blmn[0][0].size(), std::complex<double>{})));
    const auto dx = a / N_grid;
    const auto dy = b / N_grid;
    const auto dz = c / N_grid;
    const auto dV = dx * dy * dz;
    const auto N = static_cast<int>(N_grid * N_grid / size);
    for(int i = rank * N; i < (rank + 1) * N; ++i)
    {
      const auto x = dx * (i / N_grid + 0.5);
      const auto y = dy * (i % N_grid + 0.5);
      const auto xsq = x * x;
      const auto ysq = y * y;
      const auto phi = atan2(y, x);
      const double phi_quad[] = {phi, - phi, Pi + phi, Pi - phi};
      for(int k = 0; k < N_grid; ++k)
      {
        const auto z = dz * (k + 0.5);
        const auto zsq = z * z;
        const auto r = sqrt(xsq + ysq + zsq);
        const auto rn = r / rs;
        const auto xi = (rn - 1.) / (rn + 1.);
        const auto theta = acos(z / r);
        const double theta_sn[] = {theta, Pi - theta};
        const auto rho = densityxyz(x, y, z);
        for(int l = 2; l <= lmax; l += 2)
        {
          std::vector<double> rho_Uln_Iln(nmax + 1, 0.);
          for(int n = 0; n <= nmax; ++n)
          {
            rho_Uln_Iln[n] = rho * (- TSPi * pow(rn, l) / pow(1. + rn, l + l + 1)
              * RC::gegenbauer(n, l + l + 1.5, xi)) * Iln_i[l][n];
          }
          for(int m = - l; m <= l; m += 2)
          {
            if(m != 0)
            {
              for(auto theta_ : theta_sn)
              {
                for(auto phi_ : phi_quad)
                {
                  const auto Ylm_conj = std::conj(boost::math::spherical_harmonic(l, m, theta_, phi_));
                  for(int n = 0; n <= nmax; ++n)
                  {
                    blmn[l][m+lmax][n] += Ylm_conj * rho_Uln_Iln[n];
                  }
                }
              }
            }
          }
        }
      }
    }

    for(int l = 2; l <= lmax; l += 2)
    {
      for(int m = - l; m <= l; m += 2)
      {
        if(m != 0)
        {
          for(int n = 0; n <= nmax; ++n)
          {
            blmn[l][m+lmax][n] *= dV / M;
            //std::cout << "b" << l << m << n << "=" << blmn[l][m+lmax][n] << std::endl;
          }
        }
      }
    }
    return blmn;
  }
}

