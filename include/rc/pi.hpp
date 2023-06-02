/**
 * \file pi.hpp
 * \brief Contains value of Pi and various related values.
 * \author Walter Dehnen
 * 
 * C++ code written by Walter Dehnen, 1994-2005
 * Taken from falcON package and supplemented by Rimpei Chiba
 * Oxford University, Department of Physics, Theoretical Physics.
 */

#ifndef _PI_H_
#define _PI_H_

const double Pi   = 3.14159265358979323846264338328;  // Pi
const double Pih  = 0.5  * Pi;                        // Pi/2
const double Piq  = 0.25 * Pi;                        // Pi/4
const double Pi3h = 3.   * Pih;                       // 3*Pi/2
const double TPi  = 2.   * Pi;                        // 2*Pi
const double FPi  = 4.   * Pi;                        // 4*Pi
const double FPit = 4.   * Pi/3.;                     // 4*Pi/3
const double iPi  = 1.   / Pi;                        // 1/Pi
const double iPih = 1.   / Pih;                       // 2/Pi
const double iTPi = 1.   / TPi;                       // 1/(2*Pi)
const double iFPi = 1.   / FPi;                       // 1/(4*Pi)
const double iFPit= 1.   / FPit;                      // 3/(4*Pi)
const double sqPi = Pi * Pi;                          // Pi^2
const double SPi  = 1.772453850905516027298167483341; // Sqrt[Pi]
const double TSPi = 2.   * SPi;                       // 2*Sqrt[Pi]
const double STPi = 2.506628274631000502415765284811; // Sqrt[2 Pi]
const double iSPi = 1.   / SPi;                       // 1/Sqrt[Pi]

#endif
