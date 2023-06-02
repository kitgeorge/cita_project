/**
 * \file		shared/inc/export.hpp
 * \brief		Export data file.
 * \author	Rimpei Chiba
 * \date		2020
 */

#ifndef _EXPORT_H_
#define _EXPORT_H_

#include <stdlib.h>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <limits>
#include <math.h>
#include <ctype.h>
#include <assert.h>

#include "utils.hpp"

//
// Export data to filepath.
//
template <typename T>
void exportMatrix(const std::string &filepath, std::vector<std::vector<T>> &matrix);

// y over x
template <typename T>
void exportxy(const std::string &filepath, std::vector<T> &vec1D, 
double xmin, double dx, double xunit, double yunit);

// z over (x,y)
template <typename T>
void exportxyz(const std::string &filepath, std::vector<std::vector<T>> &vec2D, 
double xmin, double ymin, double dx, double dy, 
double xunit, double yunit, double zunit);

// v over (x,y,z)
template <typename T>
void exportxyzv(const std::string &filepath, std::vector<std::vector<std::vector<T>>> &vec3D, 
double xmin, double ymin, double zmin, double dx, double dy, double dz, 
double xunit, double yunit, double zunit, double vunit);

#endif