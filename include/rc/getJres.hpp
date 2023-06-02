/**
 * \file		shared/inc/getJres.hpp
 * \brief		Get Jsres.
 * \author	Rimpei Chiba
 * \date		2020
 */

#ifndef _GETJRES_H_
#define _GETJRES_H_

#pragma once
#include <stdlib.h>
#include <cstdio>
#include <vector>
#include <algorithm>
#include <array>
#include <boost/math/differentiation/autodiff.hpp>

#include "utils.hpp"
#include "units.hpp"
#include "mapXVtoAA2D.hpp"

double getJsres(int Nr, int Npsi, double wp, double Jf1, double Jf2, 
  double Jr_max, double L_max, double dJr, double dL, 
  std::vector<std::vector<double>> &wp_J, bool &gotvalue);

double getJsres(int Nr, int Npsi, double wp, double Jf1, double Jf2, 
  double Jr_max, double L_max, double dJr, double dL, 
  std::vector<std::vector<double>> &wp_J);

double getJrres(double wp, double Lres, double Jr_max, double dJr, double dL, 
  std::vector<std::vector<double>> &wp_J, bool &gotvalue);

double getJrres(double wp, double Lres, double Jr_max, double dJr, double dL, 
  std::vector<std::vector<double>> &wp_J);

#endif

