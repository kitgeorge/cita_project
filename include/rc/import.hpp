/**
 * \file		shared/inc/import.hpp
 * \brief		Import data file.
 * \author	Rimpei Chiba
 * \date		2020
 */

#ifndef _IMPORT_H_
#define _IMPORT_H_

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

// Import data stored in filepath.
//
// # comments 
// z01 z02 z03 ... 
// z11 z12 z13 ... 
// z21 z22 z23 ... 
// ...
//

void import(const std::string &filepath, std::vector<std::vector<double>> &matrix);
void importCSV(const std::string &filepath, std::vector<std::vector<double>> &matrix);

//
// struct Importxy
// Import 2D data y(x) stored in filepath.
// The data must be stored in the following format.
// The interval must be constant (i.e. dx = xi+1 - xi = const).
//
// # comments 
// x0 y0 
// x1 y1
// x2 y2
// ...
//

struct Importxy
{
  public:
    int Nx;
		double xmin, xmax, dx;

    template <typename T>
    double gety(double x, const std::vector<T> &xy) const;
    template <typename T>
		Importxy(const std::string &filepath, std::vector<T> &xy, 
      double xunit, double yunit);
};

//
// struct Importxyz
// Import 3D data z(x,y) stored in filepath.
// The data must be stored in the following format.
// The interval must be constant (i.e. dx = xi+1 - xi = const).
//
// # comments 
// x0 y0 z00
// x0 y1 z01
// x0 y2 z02
// ...
// 
// x1 y0 z10
// x1 y1 z11
// x1 y2 z12
// ...
//
// ... ... ...
//

struct Importxyz
{
  public:
    int Nx, Ny;
		double xmin, xmax, dx, ymin, ymax, dy;

    template <typename T>
    double getz(double x, double y, const std::vector<std::vector<T>> &xyz) const;
    template <typename T>
		Importxyz(const std::string &filepath, std::vector<std::vector<T>> &xyz, 
      double xunit, double yunit, double zunit);
};

#endif