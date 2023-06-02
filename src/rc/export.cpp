/**
 * \file		shared/inc/export.cpp
 * \brief		Esport data file.
 * \author	Rimpei Chiba
 * \date		2020
 */

#include "export.hpp"

template <typename T>
void exportMatrix(const std::string &filepath, std::vector<std::vector<T>> &matrix)
{
  std::ofstream fOut(filepath);
  const auto N = matrix.size();
  for(size_t i = 0; i < N; ++i)
  {
    const auto M = matrix[i].size();
    for(size_t j = 0; j < M; ++j)
    {
      fOut << matrix[i][j] << " ";
    }
    fOut << std::endl;
  }
  fOut.close();
}
template void exportMatrix(const std::string &filepath, std::vector<std::vector<int>> &matrix);
template void exportMatrix(const std::string &filepath, std::vector<std::vector<double>> &matrix);

template <typename T>
void exportxy(const std::string &filepath, std::vector<T> &vec1D, 
double xmin, double dx, double xunit, double yunit)
{
  std::ofstream fOut(filepath);
  const int Nx = vec1D.size();
  for(int i = 0; i < Nx; ++i)
  {
    const double x = (xmin + dx * i) * xunit;
    const double y = vec1D[i] * yunit;
    fOut << x << " " << y << std::endl;
  }
  fOut.close();
}
template void exportxy(const std::string &filepath, std::vector<double> &vec1D, 
double xmin, double dx, double xunit, double yunit);
template void exportxy(const std::string &filepath, std::vector<int> &vec1D, 
double xmin, double dx, double xunit, double yunit);

template <typename T>
void exportxyz(const std::string &filepath, std::vector<std::vector<T>> &vec2D, 
double xmin, double ymin, double dx, double dy, 
double xunit, double yunit, double zunit)
{
  std::ofstream fOut(filepath);
  const int Nx = vec2D.size();
  const int Ny = vec2D[0].size();
  for(int i = 0; i < Nx; ++i)
  {
    const double x = (xmin + dx * i) * xunit;
    for(int j = 0; j < Ny; ++j)
    {
      const double y = (ymin + dy * j) * yunit;
      const double z = vec2D[i][j] * zunit;
      fOut << x << " " << y << " " << z << std::endl;
    }
    fOut << std::endl;
  }
  fOut.close();
}
template void exportxyz(const std::string &filepath, std::vector<std::vector<double>> &vec2D, 
double xmin, double ymin, double dx, double dy, double xunit, double yunit, double zunit);
template void exportxyz(const std::string &filepath, std::vector<std::vector<int>> &vec2D, 
double xmin, double ymin, double dx, double dy, double xunit, double yunit, double zunit);


template <typename T>
void exportxyzv(const std::string &filepath, std::vector<std::vector<std::vector<T>>> &vec3D, 
double xmin, double ymin, double zmin, double dx, double dy, double dz, 
double xunit, double yunit, double zunit, double vunit)
{
  std::ofstream fOut(filepath);
  const int Nx = vec3D.size();
  const int Ny = vec3D[0].size();
  const int Nz = vec3D[0][0].size();
  for(int i = 0; i < Nx; ++i)
  {
    const double x = (xmin + dx * i) * xunit;
    for(int j = 0; j < Ny; ++j)
    {
      const double y = (ymin + dy * j) * yunit;
      for(int k = 0; k < Nz; ++k)
      {
        const double z = (zmin + dz * k) * zunit;
        const double v = vec3D[i][j][k] * vunit;
        fOut << x << " " << y << " " << z << " " << v << std::endl;
      }
      fOut << std::endl;
    }
  }
  fOut.close();
}
template void exportxyzv(const std::string &filepath, std::vector<std::vector<std::vector<double>>> &vec3D, 
double xmin, double ymin, double zmin, double dx, double dy, double dz, 
double xunit, double yunit, double zunit, double vunit);
template void exportxyzv(const std::string &filepath, std::vector<std::vector<std::vector<int>>> &vec3D, 
double xmin, double ymin, double zmin, double dx, double dy, double dz, 
double xunit, double yunit, double zunit, double vunit);


