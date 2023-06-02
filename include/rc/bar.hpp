/**
 * \file		shared/inc/bar.hpp
 * \brief		Functions for bar models.
 * \author	Rimpei Chiba
 * \date		2020-
 */

#ifndef _BAR_H_
#define _BAR_H_

#include <stdlib.h>
#include <math.h>
#include <boost/math/special_functions/spherical_harmonic.hpp>

#include "galaxy.hpp"
#include "force.hpp"
#include "utils.hpp"
#include "math_RC.hpp"
#include "coordinate.hpp"

namespace RC
{
  double BarGrowthDehnen(double t, double Tg); // Dehnen 1999
  double BarGrowthCosine(double t, double Tg);

  class Bar
  {    
    private:
      double Tg, Tt, wb0, eta, t2, wb0_eta, C_tr_sq, C_tr, C_tr_i;

    public:
      double wb{}, rCR{}, phib{}, phib0{}, Ib{}, Mb{};
      std::vector<std::vector<std::vector<std::complex<double>>>> blmn{};
      virtual double getwb(double t);
      virtual double getdwbdt(double t);
      virtual double getphib(double t);
      virtual void   setphib0(double phib, double t);
      virtual double Phim(double r);
      virtual double dPhimdr(double r);
      virtual double Philm(double r);
      virtual double density(const Polar &pos);
      virtual double pot(const Polar &pos);
      virtual RC::Force force(const Polar &pos);
      virtual void update(double new_wb, double new_rCR);
      virtual std::vector<std::vector<std::vector<std::complex<double>>>>
        calcblmn(int rank, int size);
      void setEvolPara(double Tg, double Tt, double wb0, double eta);
  };
  
  // Spherical harmonics Y22
  class Y22Bar : public Bar
  {
    private:
      int m;
      double b, Phim_C;
    
    public:
      double Phim(double r);
      double dPhimdr(double r);
      double Philm(double r);
      double density(const Polar &pos);
      double pot(const Polar &pos);
      RC::Force force(const Polar &pos);
      Y22Bar(int m, double A, double b);
  };

  // Ferrers 1877
  class FerrersBar : public Bar
  {
    private:
      // Bar parameters that are kept constant
      const int p;
      const double rCR_a, b_a, c_a, rhob0, C_Ib;
      // Bar parameters that are updated
      double a, b, c;
      // Parameters for basis functions
      const int lmax, nmax, N_grid;
      const double rs, M, vgsq;
      // I_ln for Hernquist & Ostriker basis function 
      std::vector<std::vector<double>> Iln_i;
      double densityxyz(double x, double y, double z);
    
    public:
      double density(const Polar &pos);
      double pot(const Polar &pos);
      RC::Force force(const Polar &pos);
      void update(double new_wb, double new_rCR);
      std::vector<std::vector<std::vector<std::complex<double>>>>
        calcblmn(int rank, int size);
      FerrersBar(int p, double a, double b, double c, double Mb, double rCR,
        int lmax, int nmax, double rs, double M);
  };
}

#endif