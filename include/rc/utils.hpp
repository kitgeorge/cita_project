/**
 * \file		shared/inc/utils.hpp
 * \brief		Some useful functions.
 * \author	Rimpei Chiba
 * \date		2018-
 */

#ifndef _UTILS_H_
#define _UTILS_H_

#include <stdlib.h>
#include <cstdlib>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <numeric>
#include <cassert>
#include <cmath>
#include <ctime>
#include <limits>
#include <sys/stat.h>
#include <stdexcept>

#include "pi.hpp"

// Output 'check'. For debugging.
inline void check()
{
  std::cout << "check" << std::endl;
}

// Output value. For debugging.
template <typename T>
inline void show(const T &value)
{
  std::cout << value << std::endl;
}

// Output name = value. For debugging.
template <typename T>
inline void show(const std::string &name, const T &value)
{
  std::cout << name << " = " << value << std::endl;
}

inline void outputCalcStartDate(int argc, char **argv, std::time_t start_date)
{
  std::cout << "\nCalculation start date : " << std::ctime(&start_date);
  for(int i = 0; i < argc; ++i) std::cout << argv[i] << " ";
  std::cout << std::endl;
}

inline void writeLogEasy(int argc, char **argv, std::string out_path)
{
  std::ofstream log(out_path + "log.txt");
  for(int i = 0; i < argc; ++i) log << argv[i] << " ";
  log << "\n\n";
  log.close();
}

inline void writeLog(int argc, char **argv, std::clock_t start_time, 
std::clock_t end_time, std::time_t start_date, std::time_t end_date, std::string out_path)
{
  std::ofstream log(out_path + "log.txt");
  std::string exec = argv[0];
  std::size_t dotslash = exec.find("./");
  exec = exec.substr(dotslash+2);
  log << "./" << exec << " ";
  for(int i = 1; i < argc; ++i) log << argv[i] << " ";
  log << "\n\n";

  double t = (double)(end_time - start_time) / CLOCKS_PER_SEC;
  double t_h = t / 3600;
  log << "Calculation start date : " << std::ctime(&start_date);
  log << "Calculation end date   : " << std::ctime(&end_date);
  log << "Total calculation time : " << t << " sec = " << t_h << " hour" << std::endl;
  log.close();

  std::cout << "Total calculation time = " << t << " sec = " << t_h << " hour\n\n";
}

// Check whether a file/directory in 'path' exsit or not
// returns 1(true) if file exist and returns 0(false) if not
inline bool fileExists(const std::string& path)
{
  // struct stat
  struct stat buffer;
  // syntax of function 'stat'
  // int stat(const char *path, struct stat *buffer)
  // path: Pointer to a string containing the path of existing file or directory.
  // buffer: Pointer to structure that stores results.
  // Returns 0 if the file-status information is obtained (file exists).
  // Returns -1 indicates that the filename or path could not be found.
  return (stat(path.c_str(), &buffer) == 0);
  // This function returns 1(true) if file exists. Otherwise 0.
}

inline bool directoryExists(const std::string& path)
{
  struct stat buffer;
  stat(path.c_str(), &buffer);
  // The macro S_ISDIR accepts stat.st_mode param and returns
  // a non-zero integer if given file is a directory, otherwise zero.
  if(stat(path.c_str(), &buffer) == 0 && S_ISDIR(buffer.st_mode)) return 1;
  return 0;
}

// Check whether a string is blank (empty OR just white space) or not
// returns 1(true) if string is blank and returns 0(false) if not
inline bool is_blank(const std::string &str)
{
  for(const char& c : str)
    if(!std::isspace(c)) return false;
  return true;
}

// Make directory 'Output' if not exist and 
// then make output directory 'out_dir' inside Output.
inline void makeOutputDirectory(const std::string& out_name)
{
  if(!directoryExists("Output"))
  {
    const auto i = mkdir("Output", ACCESSPERMS);
    if(i) show("Warning [shared/inc/utils.hpp]: 'mkdir Output' failed.");
  }
  std::string out_dir = "Output/" + out_name;
  const auto i = mkdir(out_dir.c_str(), ACCESSPERMS);
  if(i) return;
}

// Express angles within [0:2pi]
inline double Oto2PI(double x)
{
  if(std::isnan(x))
  {
    //show("Warning [shared/inc/utils.hpp]: Argument const double x of function Oto2PI is nan.");
    return 0;
  }
  else
  {
    // Express x in [0,2pi)
    auto a = fmod(TPi + fmod(x, TPi), TPi);
    assert((a >= 0.0) && (a <= TPi));
    return a;
  }
}

// Express angles within [-pi:pi]
inline double PItoPI(double x)
{
  // Express x in [0,2pi]
  auto a = Oto2PI(x);
  // Express x in [-pi,pi]
  if(a > Pi) a -= TPi;
  assert((- Pi <= a) && (a <= Pi));
  return a;
}

// Functions to receive command line arguments. 
// Exit the process if type is mismatched.
inline int recv_i(const std::string &argName, const std::string &str)
{
  int i;
  try
  {
    i = std::stoi(str);
  }
  catch(const std::invalid_argument& e)
  {
    std::cerr << "\nOption error: Value of '" << argName 
              << "' must be integer.\n" << std::endl;
    exit(EXIT_FAILURE);
  }
  return i;
}

inline size_t recv_zu(const std::string &argName, const std::string &str)
{
  std::istringstream iss(str);
  size_t size;
  iss >> size;
  if(iss.fail())
  {
    std::cerr << "\nOption error: Value of '" << argName 
              << "' must be integer.\n" << std::endl;
    exit(EXIT_FAILURE);
  }
  return size;
}

inline double recv_f(const std::string &argName, const std::string &value)
{
  double f;
  try
  {
    f = std::stof(value);
  }
  catch(const std::invalid_argument& e)
  {
    std::cerr << "\nOption error: Value of '" << argName
              << "' must be double.\n" << std::endl;
    exit(EXIT_FAILURE);
  }
  return f;
}

inline bool recv_b(const std::string &argName, const std::string &str)
{
  bool b;
  std::istringstream is(str);
  // first try simple integer conversion
  is >> b;
  if(is.fail())
  {
    // simple integer failed; try boolean
    is.clear();
    is >> std::boolalpha >> b;
  }
  if(is.fail())
  {
    std::cerr << "\nOption error: Value of '" << argName 
              << "' must be bool.\n" << std::endl;
    exit(EXIT_FAILURE);
  }
  return b;
}

inline std::string fill0(int i, int width)
{
  std::stringstream ss;
  ss << std::setfill('0') << std::setw(width) << i;
  return ss.str();
}

inline std::string fixedDecimal(double f, int digit)
{
  std::stringstream ss;
  ss << std::fixed << std::setprecision(digit) << f;
  return ss.str();
}

std::string C16bit(const double x);
std::string RGB(const double R, const double G, const double B);
std::string K2LB(const double x);
std::string K2LR(const double x);
double calcMean(const std::vector<double> &sample);
double calcStdv_biased(const std::vector<double> &sample, const double mean);
double calcStdv_unbiased(const std::vector<double> &sample, const double mean);


#endif


