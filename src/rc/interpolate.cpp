/**
 * \file		shared/inc/interpolate.cpp
 * \brief		Functions to interpolate data.
 * \author	Rimpei Chiba
 * \date		2020
 */

#include "interpolate.hpp"

template <typename T>
std::vector<double> linearFit(const std::vector<std::vector<T>> &xy)
{
  std::vector<T> x,y;
  for(auto &vec : xy)
  {
    x.push_back(vec[0]);
    y.push_back(vec[1]);
  }
  const auto n    = xy.size();
  const auto s_x  = std::accumulate(x.begin(), x.end(), 0.);
  const auto s_y  = std::accumulate(y.begin(), y.end(), 0.);
  const auto s_xx = std::inner_product(x.begin(), x.end(), x.begin(), 0.);
  const auto s_xy = std::inner_product(x.begin(), x.end(), y.begin(), 0.);
  const auto a    = (n * s_xy - s_x * s_y) / (n * s_xx - s_x * s_x);
  const auto b    = (s_y - a * s_x) / n;
  std::vector<double> res;
  res.push_back(a);
  res.push_back(b);
  return res;
}
template std::vector<double> linearFit(const std::vector<std::vector<int>>& xy);
template std::vector<double> linearFit(const std::vector<std::vector<float>>& xy);
template std::vector<double> linearFit(const std::vector<std::vector<double>>& xy);

void bilinearInterpolation(const std::string &file_path, int scale_x, int scale_y)
{
  std::vector<std::vector<double>> f;
  Importxyz in(file_path,f,1,1,1);
  if(!f.empty())
  {
    const auto dx_i = 1. / in.dx;
    const auto dy_i = 1. / in.dy;
    const auto dxy_i = dx_i * dy_i;
    const auto Nx_m1 = in.Nx - 1; 
    const auto Ny_m1 = in.Ny - 1;
    const auto dx_new = in.dx / scale_x;
    const auto dy_new = in.dy / scale_y;

    // Create new file with name "<original name>_bi.dat"
    std::size_t ext = file_path.find(".dat");
    std::string file_path_new = file_path.substr(0, ext) + "_bi.dat";
    std::ofstream fOut(file_path_new);
    for(int i = 0; i < Nx_m1; ++i)
    {
      const auto x1 = in.xmin + in.dx * i;
      const auto x2 = x1 + in.dx;
      for(int s = 0; s < scale_x; ++s)
      {
        const auto x = x1 + dx_new * s;
        const auto x2_x = x2 - x;
        const auto x_x1 = x - x1;
        for(int j = 0; j < Ny_m1; ++j)
        {
          const auto y1 = in.ymin + in.dy * j;
          const auto y2 = y1 + in.dy;
          for(int t = 0; t < scale_y; ++t)
          {
            const auto y = y1 + dy_new * t;
            const auto y2_y = y2 - y;
            const auto y_y1 = y - y1;
            const auto f_bi = dxy_i*(y2_y*(x2_x*f[i][j]+x_x1*f[i+1][j])
                              + y_y1*(x2_x*f[i][j+1]+x_x1*f[i+1][j+1]));
            fOut << x << " " << y << " " << f_bi << std::endl;
          }
        }
        const auto f_bi = dx_i*(x2_x*f[i][Ny_m1]+x_x1*f[i+1][Ny_m1]);
        fOut << x << " " << in.ymax << " " << f_bi << std::endl;
        fOut << std::endl;
      }
    }
    for(int j = 0; j < Ny_m1; ++j)
    {
      const auto y1 = in.ymin + in.dy * j;
      const auto y2 = y1 + in.dy;
      for(int t = 0; t < scale_y; ++t)
      {
        const auto y = y1 + dy_new * t;
        const auto y2_y = y2 - y;
        const auto y_y1 = y - y1;
        const auto f_bi = dy_i*(y2_y*f[Nx_m1][j]+y_y1*f[Nx_m1][j+1]);
        fOut << in.xmax << " " << y << " " << f_bi << std::endl;
      }
    }
    fOut << in.xmax << " " << in.ymax << " " 
         << f[Nx_m1][Ny_m1] << std::endl;
    fOut.close();
  }
}