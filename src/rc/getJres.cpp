/**
 * \file		shared/src/getJres.cpp
 * \brief		Get Jsres.
 * \author	Rimpei Chiba
 * \date		2020
 */

#include "getJres.hpp"

//
// Function to get Jsres from (wp,Jf1,Jf2)
//

namespace RC {

double getJsres(int Nr, int Npsi, double wp, double Jf1, double Jf2, 
  double Jr_max, double L_max, double dJr, double dL, 
  std::vector<std::vector<double>> &wp_J, bool &gotvalue)
{
  const bool report = false;
  const auto dJs = L_max * 0.001;
  double Js_min{}, Js_max{};
  if(Nr > 0)
  {
    if(Npsi > 0)
    {
      Js_min = std::max(- Jf1 / Nr, - Jf2 / Npsi);            
      Js_max = std::min((Jr_max - Jf1) / Nr, (L_max - Jf2) / Npsi);
    }
    else if(Npsi < 0)
    {
      Js_min = std::max(- Jf1 / Nr, (L_max - Jf2) / Npsi);
      Js_max = std::min((Jr_max - Jf1) / Nr, - Jf2 / Npsi);
    }
    else // if Npsi == 0
    {
      Js_min = - Jf1 / Nr;            
      Js_max = (Jr_max - Jf1) / Nr;
    }
  }
  else if(Nr < 0)
  {
    if(Npsi > 0)
    {
      Js_min = std::max((Jr_max - Jf1) / Nr, - Jf2 / Npsi);            
      Js_max = std::min(- Jf1 / Nr, (L_max - Jf2) / Npsi);
    }
    else if(Npsi < 0)
    {
      Js_min = std::max((Jr_max - Jf1) / Nr, (L_max - Jf2) / Npsi);            
      Js_max = std::min(- Jf1 / Nr, - Jf2 / Npsi);
    }
    else // if Npsi == 0
    {
      Js_min = (Jr_max - Jf1) / Nr;
      Js_max = - Jf1 / Nr;            
    }
  }
  else // if Nr == 0
  {
    if(Npsi > 0)
    {
      Js_min = - Jf2 / Npsi;
      Js_max = (L_max - Jf2) / Npsi;
    }
    else if(Npsi < 0)
    {
      Js_min = (L_max - Jf2) / Npsi;
      Js_max = - Jf2 / Npsi;
    }
  }
  Js_min += dJs;
  Js_max -= dJs;

  // Find Js=Jsres where wp(Jsres)-wp = 0 using bisection method
  const auto dwp_err = 1.e-4 * Units::Gyr_i; 
  const int count_max = 100;
  double Js0 = Js_min;
  const auto Jr0 = Jf1 + Nr * Js0;
  const auto L0  = Jf2 + Npsi * Js0;
  const auto wp0 = RC::XGivenJ(Jr0, L0, dJr, dL, wp_J);
  double dwp0 = wp0 - wp;
  double Js1 = Js_max;
  const auto Jr1 = Jf1 + Nr * Js1;
  const auto L1  = Jf2 + Npsi * Js1;
  const auto wp1 = XGivenJ(Jr1, L1, dJr, dL, wp_J);
  double dwp1 = wp1 - wp;
  if(report)
  {
    show("Jsmin", Js_min * Units::kpc2Gyr_i);
    show("Jsmax", Js_max * Units::kpc2Gyr_i);
    show("Jr0", Jr0 * Units::kpc2Gyr_i);
    show("L0 ", L0 * Units::kpc2Gyr_i);
    show("Jr1", Jr1 * Units::kpc2Gyr_i);
    show("L1 ", L1 * Units::kpc2Gyr_i);
    show("wp ", wp  * Units::Gyr);
    show("wp0", wp0 * Units::Gyr);
    show("wp1", wp1 * Units::Gyr);
    show("dwp0=wp0-wp", dwp0 * Units::Gyr);
    show("dwp1=wp1-wp", dwp1 * Units::Gyr);
    std::cout << std::endl;
  }

  if(dwp0 * dwp1 < 0)
  {
    double dwp = dwp0;
    double Js{};
    int count{};
    while(fabs(dwp) > dwp_err)
    {
      Js = (Js0 + Js1) * 0.5;
      const auto Jr = Jf1 + Nr * Js;
      const auto L  = Jf2 + Npsi * Js;
      dwp = XGivenJ(Jr, L, dJr, dL, wp_J) - wp;
      if(dwp0 * dwp > 0)
      {
        Js0 = Js;
        dwp0 = dwp;
      }
      else if(dwp * dwp1 > 0)
      {
        Js1 = Js;
        dwp1 = dwp;
      }
      else break; // dwp == 0
      ++ count;
      if(count > count_max)
      {
        show("Warning [shared/src/getJsres.cpp]: count > count_max at (dwp0 * dwp1 < 0)");
        show("Jf1", Jf1 * Units::kpc2Gyr_i);
        show("Jf2", Jf2 * Units::kpc2Gyr_i);
        show("wp", (dwp+wp) * Units::Gyr);
        gotvalue = false;
        return 0.;
      }
      if(report) show("wp", (dwp+wp) * Units::Gyr);
    }
    gotvalue = true;
    return Js;
  }
  else // happens at large (Jr,L)
  {
    if(wp0 < wp)
    {
      const auto dJs = 1. * Units::kpckms;
      const auto Js_cut = 200. * Units::kpckms;
      double Js = Js_min;
      double dwp{};
      bool cross_res = false;
      while(Js < Js_cut)
      {
        const auto Jr = Jf1 + Nr * Js;
        const auto L  = Jf2 + Npsi * Js;
        if((0 <= Jr) && (Jr < Jr_max) 
        && (0 <= L) && (L < L_max))
        {
          dwp = XGivenJ(Jr, L, dJr, dL, wp_J) - wp;
          if(dwp > 0)
          {
            cross_res = true;
            break;
          }
          Js += dJs;
        }
        else
        {
          gotvalue = false;
          return 0.;
        }
      }
      if(cross_res)
      {
        Js0 = Js;
        dwp0 = dwp;
        int count{};
        while(fabs(dwp) > dwp_err)
        {
          Js = (Js0 + Js1) * 0.5;
          const auto Jr = Jf1 + Nr * Js;
          const auto L  = Jf2 + Npsi * Js;
          dwp = XGivenJ(Jr, L, dJr, dL, wp_J) - wp;
          if(dwp0 * dwp > 0)
          {
            Js0 = Js;
            dwp0 = dwp;
          }
          else if(dwp * dwp1 > 0)
          {
            Js1 = Js;
            dwp1 = dwp;
          }
          else break; // dwp == 0
          ++ count;
          if(count > count_max)
          {
            show("Warning [shared/src/getJsres.cpp]: count > count_max at (dwp0 * dwp1 > 0)");
            show("Jf1", Jf1 * Units::kpc2Gyr_i);
            show("Jf2", Jf2 * Units::kpc2Gyr_i);
            show("wp",(dwp+wp) * Units::Gyr);
            gotvalue = false;
            return 0.;
          }
          if(report) show("wp",(dwp+wp)*Units::Gyr);
        }
        gotvalue = true;
        return Js;
      }
      else
      {
        gotvalue = false;
        return 0;
      }
    }
    else
    {
      gotvalue = false;
      return 0;
    }
  }
}

double getJsres(int Nr, int Npsi, double wp, double Jf1, double Jf2, 
  double Jr_max, double L_max, double dJr, double dL, 
  std::vector<std::vector<double>> &wp_J)
{
  bool gotvalue;
  return getJsres(Nr, Npsi, wp, Jf1, Jf2, Jr_max, L_max, dJr, dL, wp_J, gotvalue);
}

double getJrres(double wp, double Lres, double Jr_max, 
  double dJr, double dL, std::vector<std::vector<double>> &wp_J, bool &gotvalue)
{
  // Find Jr=Jrres where wp(Lres,Jrres)-wp = 0 using bisection method
  if(Lres < 0)
  {
    gotvalue = false;
    return 0.;
  }
  else
  {
    const bool report = false;
    const auto dwp_err = 1.e-4 * Units::Gyr_i; 
    const int count_max = 100;
    double Jr0 = 0.;
    const auto wp0 = RC::XGivenJ(Jr0, Lres, dJr, dL, wp_J);
    auto dwp0 = wp0 - wp;
    auto Jr1 = Jr_max;
    const auto wp1 = XGivenJ(Jr1, Lres, dJr, dL, wp_J);
    auto dwp1 = wp1 - wp;
    if(report)
    {
      show("wp0", wp0 * Units::Gyr);
      show("wp1", wp1 * Units::Gyr);
    }

    if(dwp0 * dwp1 > 0)
    {
      gotvalue = false;
      return 0.;
    }
    else
    {
      auto dwp = dwp0;
      double Jr{};
      int count{};
      while(fabs(dwp) > dwp_err)
      {
        Jr = (Jr0 + Jr1) * 0.5;
        dwp = XGivenJ(Jr, Lres, dJr, dL, wp_J) - wp;
        if(dwp0 * dwp > 0)
        {
          Jr0 = Jr;
          dwp0 = dwp;
        }
        else if(dwp * dwp1 > 0)
        {
          Jr1 = Jr;
          dwp1 = dwp;
        }
        else break; // dwp == 0
        ++ count;
        if(count > count_max)
        {
          show("Warning [shared/src/getJres.cpp]: count > count_max at (dwp0 * dwp1 < 0)");
          show("Jr", Jr * Units::kpc2Gyr_i);
          show("Lres", Lres * Units::kpc2Gyr_i);
          show("wp", (dwp + wp) * Units::Gyr);
          gotvalue = false;
          return 0.;
        }
        if(report) show("wp", (dwp + wp) * Units::Gyr);
      }
      gotvalue = true;
      return Jr;
    }
  }
}

double getJrres(double wp, double Lres, double Jr_max, double dJr, double dL, 
  std::vector<std::vector<double>> &wp_J)
{
  bool gotvalue;
  return getJrres(wp, Lres, Jr_max, dJr, dL, wp_J, gotvalue);
}

}







