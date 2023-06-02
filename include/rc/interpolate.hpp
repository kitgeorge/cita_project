/**
 * \file		shared/inc/interpolate.hpp
 * \brief		Functions to interpolate data.
 * \author	Rimpei Chiba
 * \date		2020
 */

#ifndef _INTERPOLATE_H_
#define _INTERPOLATE_H_

#include <stdlib.h>
#include <cstdlib>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <numeric>

#include "import.hpp"

// Linear regression
template <typename T>
std::vector<double> linearFit(const std::vector<std::vector<T>> &xy);

// Generate fine data using bilinear interpolation
// (for plotting contour maps in gnuplot). 
void bilinearInterpolation(const std::string &file_path, int scale_x, int scale_y);

#endif