/**
 * \file		shared/src/utils.cpp
 * \brief		Some useful functions. 
 * \author	Rimpei Chiba
 * \date		2020 
 */

#include "utils.hpp"

std::string C16bit(const double x)
{
  const auto a = lrint(x * 255);
  const int x16_1 = a / 16;
  std::string c1 =
    x16_1 < 10 ? std::to_string(x16_1) : 
    x16_1 == 10 ? "A" : 
    x16_1 == 11 ? "B" : 
    x16_1 == 12 ? "C" : 
    x16_1 == 13 ? "D" : 
    x16_1 == 14 ? "E" : "F";
  const int x16_2 = a % 16;
  std::string c2 =
    x16_2 < 10 ? std::to_string(x16_2) : 
    x16_2 == 10 ? "A" : 
    x16_2 == 11 ? "B" : 
    x16_2 == 12 ? "C" : 
    x16_2 == 13 ? "D" : 
    x16_2 == 14 ? "E" : "F";
  return c1 + c2;
}

std::string RGB(const double R, const double G, const double B)
{
  return C16bit(R) + C16bit(G) + C16bit(B);
}

std::string K2LB(const double x)
{
  const auto R = x < 0.5 ? 0 : 0.7 * (x - 0.5) * 2;
  const auto G = x < 0.5 ? 0 : 0.95 * (x - 0.5) * 2;
  const auto B = x < 0.5 ? x + x : 1.;
  return RGB(R, G, B);
}

std::string K2LR(const double x)
{
  const auto R = x < 0.5 ? x + x : 1.;
  const auto G = x < 0.5 ? 0 : 0.7 * (x - 0.5) * 2;
  const auto B = x < 0.5 ? 0 : 0.8 * (x - 0.5) * 2;
  return RGB(R, G, B);
}

double calcMean(const std::vector<double> &sample)
{
  // Designed to output nan when sample.size() == 0
  if(sample.size() == 0) return std::numeric_limits<double>::quiet_NaN();
  else return std::accumulate(sample.begin(), sample.end(), 0.) / sample.size();
}

double calcStdv_biased(const std::vector<double> &sample, const double mean)
{
  // Designed to output nan when d.o.f < 1
  if(sample.size() < 1) return std::numeric_limits<double>::quiet_NaN();
  else
  {
    double std = 0.;
    for(size_t i = 0; i < sample.size(); ++i)
    {
      auto dx = sample[i] - mean;
      std += dx * dx;
    }
    return sqrt(std / sample.size());
  }
}

double calcStdv_unbiased(const std::vector<double> &sample, const double mean)
{
  // Designed to output nan when d.o.f < 1
  if(sample.size() < 2) return std::numeric_limits<double>::quiet_NaN();
  else
  {
    double std = 0.;
    for(size_t i = 0; i < sample.size(); ++i)
    {
      auto dx = sample[i] - mean;
      std += dx * dx;
    }
    return sqrt(std / (sample.size() - 1));
  }
}

