/**
 * \file		shared/inc/import.cpp
 * \brief		Import data file and store in vector.
 * \author	Rimpei Chiba
 * \date		2020
 */

#include "import.hpp"

void import(const std::string &filepath, std::vector<std::vector<double>> &matrix)
{
  std::ifstream data(filepath);
  if(!data.is_open())
  {
    std::cout << "Warning [shared/src/import.cpp]: file \'" 
    << filepath << "\' does not exist." << std::endl;
  }
  matrix.clear();
  std::string line;
  while(std::getline(data, line))
  {
    // If line is not empty nor just filled with spaces
    // and if # is not the 0th character of line,
    // then read data
    if(!is_blank(line) && (line.find('#') != 0))
    {
      std::stringstream ss(line);
      double value;
      std::vector<double> vec;
      while(ss >> value)
      {
        vec.push_back(value);
      }
      matrix.push_back(vec);
    }
  }
}

void importCSV(const std::string &filepath, std::vector<std::vector<double>> &matrix)
{
  std::ifstream data(filepath);
  if(!data.is_open())
  {
    std::cout << "Warning [shared/src/import.cpp]: file \'" 
    << filepath << "\' does not exist." << std::endl;
  }
  matrix.clear();
  std::string line;
  while(std::getline(data, line))
  {
    if(!is_blank(line) && (line.find('#') != 0))
    {
      std::stringstream ss(line);
      std::string component;
      std::vector<double> vec;
      while(std::getline(ss,component,','))
      {
        if(component.size())
          vec.push_back(std::stof(component));
        else
          vec.push_back(0);
      }
      matrix.push_back(vec);
    }
  }
}


//
// Importxy
//

template <typename T>
Importxy::Importxy(const std::string &filepath, std::vector<T> &xy, 
double xunit, double yunit) : Nx(0)
{
  if(!fileExists(filepath))
  {
    std::cout << "Warning [shared/src/import.cpp]: file \'" 
    << filepath << "\' does not exist." << std::endl;
  }
  else
  {
    std::ifstream data(filepath);
    std::string line;

    // Find Nx
    while(std::getline(data, line))
    {
      if(!is_blank(line) && (line.find('#') != 0)) ++ Nx;
    }
    if(Nx < 2)
    {
      std::cout << "Warning [shared/src/import.cpp]: Nx of \'" 
      << filepath << "\' is " << Nx << std::endl;
    }
    data.clear();
    data.seekg(0);

    // Import y(x) into xy
    xy.resize(Nx, 0);
    for(int i = 0; i < Nx; ++i)
    {
      // Skip comments starting from #
      while(data.peek() == '#')
      {
        data.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
      }
      double x,y;
      data >> x >> y;
      if(i == 0)
      {
        // Initialize xmin,xmax
        xmin = x;
        xmax = x;
      }
      else
      {
        // Update xmin,xmax
        xmin = x < xmin ? x : xmin;
        xmax = x > xmax ? x : xmax;
      }
      // Get y(x)
      xy[i] = y * yunit;
    }
    data.close();
  }
  xmin *= xunit;
  xmax *= xunit;
  // Calculate dx
  dx = (xmax - xmin) / (Nx - 1);
}

template Importxy::Importxy(const std::string &filepath, std::vector<double> &xy, 
double xunit, double yunit);
template Importxy::Importxy(const std::string &filepath, std::vector<int> &xy, 
double xunit, double yunit);

template <typename T>
double Importxy::gety(double x, const std::vector<T> &xy) const
{
  if((xmin <= x) && (x < xmax))
  {
    const auto i = floor((x - xmin) / dx);
    const auto x1 = xmin + dx * i;
    const auto x2 = xmin + dx * (i+1);
    const T y1 = xy[i];
    const T y2 = xy[i+1];
    const double y = (y1*(x2-x) + y2*(x-x1)) / dx;
    return y;
  }
  else
  {
   std::cout << "Warning [shared/src/import.cpp]: Importxy::gety cannot get y."
   << "\nx out of range. returned y=0." << std::endl;
   return 0.;
  } 
}

template double Importxy::gety(double x, const std::vector<int> &xy) const;
template double Importxy::gety(double x, const std::vector<double> &xy) const;


//
// Importxyz
//

template <typename T>
Importxyz::Importxyz(const std::string &filepath, std::vector<std::vector<T>> &xyz, 
double xunit, double yunit, double zunit) : Nx(0), Ny(0)
{
  if(!fileExists(filepath))
  {
    std::cout << "Warning [shared/src/import.cpp]: file \'" 
    << filepath << "\' does not exist." << std::endl;
  }
  else
  {
    std::ifstream data(filepath);
    std::string line;
    // Find Ny
    while(std::getline(data, line))
    {
      if(is_blank(line) && (Ny != 0)) break;
      else if(line.find('#') != 0) ++ Ny;
    }
    data.clear();
    data.seekg(0);
    if(Ny < 2)
    {
      std::cout << "Warning [shared/src/import.cpp]: Ny of \'" 
      << filepath << "\' is " << Ny << std::endl;
    }
    // Find Nx
    while(std::getline(data, line))
    {
      if(!is_blank(line) && (line.find('#') != 0)) ++ Nx;
    }
    float nx = static_cast<float>(Nx) / Ny;
    if(abs(nx - static_cast<int>(nx)) > 0)
    {
      std::cout << "Warning [shared/src/import.cpp]: Data points z(x,y) in \'" 
      << filepath << "\' do not have a Nx*Ny grid." << std::endl;
    }
    Nx = static_cast<int>(nx);
    data.clear();
    data.seekg(0);
    if(Nx < 2)
    {
      std::cout << "Warning [shared/src/import.cpp]: Nx of \'" 
      << filepath << "\' is " << Nx << std::endl;
    }
    // Import z(x,y) into xyz
    xyz.resize(Nx, std::vector<T>(Ny,0));
    for(int i = 0; i < Nx; ++i)
    {
      for(int j = 0; j < Ny; ++j)
      {
        // Skip comments starting from #
        while(data.peek() == '#')
        {
          data.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        }
        double x,y,z;
        data >> x >> y >> z;
        if((i == 0) && (j == 0))
        {
          // Initialize xmin,xmax,ymin,ymax
          xmin = x;
          xmax = x;
          ymin = y;
          ymax = y;
        }
        else
        {
          // Update xmin,xmax,ymin,ymax
          xmin = x < xmin ? x : xmin;
          xmax = x > xmax ? x : xmax;
          ymin = y < ymin ? y : ymin;
          ymax = y > ymax ? y : ymax;
        }
        // Get z(x,y)
        xyz[i][j] = z * zunit;
      }
    }
    data.close();
    xmin *= xunit;
    xmax *= xunit;
    ymin *= yunit;
    ymax *= yunit;
    // Calculate dx,dy
    dx = (xmax - xmin) / (Nx - 1);
    dy = (ymax - ymin) / (Ny - 1);
  }
}

template Importxyz::Importxyz(const std::string &filepath, std::vector<std::vector<double>> &xyz, 
double xunit, double yunit, double zunit);
template Importxyz::Importxyz(const std::string &filepath, std::vector<std::vector<int>> &xyz, 
double xunit, double yunit, double zunit);

template <typename T>
double Importxyz::getz(double x, double y, const std::vector<std::vector<T>> &xyz) const
{
  if((xmin <= x) && (x < xmax) && (ymin <= y) && (y < ymax))
  {
    const auto i = floor((x - xmin) / dx);
    const auto j = floor((y - ymin) / dy);
    const auto x1 = xmin + dx * i;
    const auto x2 = xmin + dx * (i+1);
    const auto y1 = ymin + dy * j;
    const auto y2 = ymin + dy * (j+1);
    const T z11 = xyz[i][j];
    const T z12 = xyz[i][j+1];
    const T z21 = xyz[i+1][j];
    const T z22 = xyz[i+1][j+1];
    const double z = ((x2-x)*(y2-y)*z11 + (x-x1)*(y2-y)*z21 
                    + (x2-x)*(y-y1)*z12 + (x-x1)*(y-y1)*z22) / (dx*dy);
    return z;
  }
  else
  {
   std::cout << "Warning [shared/src/import.cpp]: Importxyz::getz cannot get z."
   << "\nxy out of range. returned z=0." << std::endl;
   return 0.;
  } 
}

template double Importxyz::getz(double x, double y, 
const std::vector<std::vector<int>> &xyz) const;
template double Importxyz::getz(double x, double y, 
const std::vector<std::vector<double>> &xyz) const;
