/**
 * \file		shared/inc/mapXVtoAA2D.hpp
 * \brief		Map (x,v) to (theta,J)
 * \author	Rimpei Chiba
 * \date		2019-
 */

// Modified to work with AxsymFuncs class

#ifndef _MAPXVTOAA2D_H_
#define _MAPXVTOAA2D_H_

#pragma once
#include <stdlib.h>
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <vector>
#include <math.h>
#include <assert.h>
#include <complex>
#include <functional>

#include "pi.hpp"
#include "units.hpp"
#include "galaxy.hpp"
#include "utils.hpp"
#include "potential.hpp"
#include "NewtonRaphson.hpp"
#include "math_RC.hpp"
#include "axsym_funcs.hpp"

namespace RC {

class MapXVtoAA2D
{
	private:
		bool warn;
		potential::AxsymFuncs *Phi;
		const double dr_apsis, dr_apsis_err, r_upper_limit, r_lower_limit, tau_max;
		void clear();
		bool findApsis(double L, double E); 

  public:
		double L, E, r_apo, r_peri, Jr, Jpsi, Tr, Tpsi, wr, wpsi, thetar, thetapsi;
		void epicycle(double L);
		bool mapLEtoJ(double L_, double E_, int N_tau);
		bool mapXVtoAA(double r, double psi, double vr, double vpsi, int N_tau);
		double W_k(int kr, int kp, const std::function<double (double)> &Philm, int N_tau);
		MapXVtoAA2D(potential::AxsymFuncs *Phi);
};

void gridEOverJ(potential::AxsymFuncs *Phi, int N_tau, int N_Jr, int N_L, 
double dJr, double dL, std::vector<std::vector<double>> &E_J, bool prog);

double XGivenJ(double Jr, double L, double dJr, double dL, std::vector<std::vector<double>> &X_J);

}

#endif

