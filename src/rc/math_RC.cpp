/**
 * \file		shared/inc/math_RC.hpp
 * \brief   Mathematic functions.
 * \author	Rimpei Chiba
 * \date		2020-
 */

#include "math_RC.hpp"

namespace RC
{
  template <typename T>
  double norm(std::vector<T> vec)
  {return sqrt(std::inner_product(vec.begin(), vec.end(), vec.begin(), 0));}
  template double norm(std::vector<int> vec);
  template double norm(std::vector<float> vec);
  template double norm(std::vector<double> vec);
  template double norm(std::vector<long double> vec);

  template <typename T>
  T sq(const T x) {return x * x;}
  template int sq(int x);
  template float sq(float x);
  template double sq(double x);
  template long double sq(long double x);

  template <typename T>
  T cb(const T x) {return x * x * x;}
  template int cb(int x);
  template float cb(float x);
  template double cb(double x);
  template long double cb(long double x);

  template <typename T>
  T sgn(const T x) {return (T(0) < x) - (x < T(0));}
  template int sgn(int x);
  template float sgn(float x);
  template double sgn(double x);
  template long double sgn(long double x);

  unsigned factorial(unsigned n)
  {
    if(n > 1) return n * factorial(n - 1);
    else return 1;
  }

  unsigned doublefactorial(unsigned n)
  {
    if(n > 1) return n * doublefactorial(n - 2);
    else return 1;
  }

  template<typename Real>
  Real gegenbauer(unsigned int n, Real a, Real x)
  {
    //
    // Gegenbauer polynomials
    //
    // C^a_n(x) = 1/n * [2*x*(n+a-1)*C^a_{n-1}(x) - (n+2*a-2)*C^a_{n-2}(x)]
    //          = 1/n * [x*(2*n+g)*C^a_{n-1}(x) - (n+g)*C^a_{n-2}(x)]
    // where  g = 2 * (a - 1)
    //
    static_assert(!std::is_integral<Real>::value, 
      "Gegenbauer polynomials required floating point arguments.");

    if(n == 0) return Real(1);
    Real g = 2 * (a - 1);
    Real C0 = 1.;
    Real C1 = 2 * a * x;
    Real Ck = C1;
    unsigned int k = 2;
    while(k <= n)
    {
      Ck = (x * (k + k + g) * C1 - (k + g) * C0) / k;
      C0 = C1;
      C1 = Ck;
      k += 1;
    }
    return Ck;
  }
  template float gegenbauer(unsigned n, float a, float x);
  template double gegenbauer(unsigned n, double a, double x);
  template long double gegenbauer(unsigned n, long double a, long double x);

  template<typename Real>
  Real gegenbauer_derivative(unsigned int n, Real a, Real x, unsigned int k)
  {
    //
    // d/dx C^a_n(x) = 2 * a * C^{a+1}_{n-1}(x)
    // (e.g. Piretzidis & Sideris 2021)
    //
    if(k > n) return Real(0);
    Real gegen = gegenbauer<Real>(n - k, a + k, x);
    Real scale = 1;
    for(unsigned j = 0; j < k; ++j)
    {
      scale *= a + a;
      a += 1;
    }
    return scale * gegen;
  }
  template float gegenbauer_derivative(unsigned int n, float a, float x, unsigned int k);
  template double gegenbauer_derivative(unsigned int n, double a, double x, unsigned int k);
  template long double gegenbauer_derivative(unsigned int n, long double a, long double x, unsigned int k);

  double WignerdMatrix(unsigned l, int mp, int m, double beta)
  {
    if((unsigned)abs(m) > l) throw std::domain_error("l >= |m| is required.");
    else if((unsigned)abs(mp) > l) throw std::domain_error("l >= |m'| is required.");
    auto cos_hb = cos(0.5 * beta);
    auto sin_hb = sin(0.5 * beta);
    auto l2 = l + l;
    auto m_mp = m - mp;

    int t_min{};
    if(m_mp > 0) t_min = 0;
    else t_min = - m_mp;

    int t_max{};
    if(m + mp < 0) t_max = l + mp;
    else t_max = l - m;

    double d{};
    for(int t = t_min; t <= t_max; ++t)
    {
      auto t2 = t + t;
      auto buf = pow(cos_hb, l2 - m_mp - t2) *
                 pow(sin_hb, t2 + m_mp) /
                 (factorial(l - m - t) *
                 factorial(l + mp - t) *
                 factorial(t) *
                 factorial(t + m_mp));
      if(t % 2 == 0) d += buf;
      else d -= buf;
    }
    d *= sqrt(factorial(l + m) *
              factorial(l - m) *
              factorial(l + mp) *
              factorial(l - mp));
    return d;
  }

  double WignerdMatrix(unsigned l, int mp, int m, double beta, double &derivative)
  {
    if((unsigned)abs(m) > l) throw std::domain_error("l >= |m| is required.");
    else if((unsigned)abs(mp) > l) throw std::domain_error("l >= |m'| is required.");
    auto cos_hb = cos(0.5 * beta);
    auto sin_hb = sin(0.5 * beta);
    auto tan_hb = sin_hb / cos_hb;
    auto cot_hb = cos_hb / sin_hb;
    auto m_mp = m - mp;

    int t_min{};
    if(m_mp > 0) t_min = 0;
    else t_min = - m_mp;

    int t_max{};
    if(m + mp < 0) t_max = l + mp;
    else t_max = l - m;

    double d{}, a{};
    for(int t = t_min; t <= t_max; ++t)
    {
      int c = 2 * (l - t) - m_mp;
      int s = 2 * t + m_mp;
      auto buf = pow(cos_hb, c) * pow(sin_hb, s) /
                 (factorial(l - m - t) * 
                  factorial(l + mp - t) *
                  factorial(t) * 
                  factorial(t + m_mp));
      auto buf_deriv = buf * (- c * tan_hb + s * cot_hb);
      if(t % 2 == 0)
      {
        d += buf;
        a += buf_deriv;
      }
      else
      {
        d -= buf;
        a -= buf_deriv;
      }
    }
    auto f = sqrt(factorial(l + m) * 
                  factorial(l - m) *
                  factorial(l + mp) * 
                  factorial(l - mp));
    d *= f;
    derivative = 0.5 * a * f;
    return d;
  }

  std::complex<double> WignerDMatrix(unsigned l, int mp, int m, 
    double alpha, double beta, double gamma)
  {
    const auto a = - mp * alpha - m * gamma;
    std::complex<double> D(cos(a), sin(a));
    return D * WignerdMatrix(l, mp, m, beta);
  }

  double WignerdMatrix22(int mp, double beta)
  {
    if((unsigned)abs(mp) > 2) throw std::domain_error("l >= |m'| is required.");
    if     (mp == 2) return pow(cos(0.5 * beta), 4);
    else if(mp == 1) return sin(beta) * pow(cos(0.5 * beta), 2);
    else if(mp == 0) return 0.6123724356957945 * pow(sin(beta), 2); // sqrt(6) / 4
    else if(mp ==-1) return sin(beta) * pow(sin(0.5 * beta), 2);
    else             return pow(sin(0.5 * beta), 4);
  }

  double WignerdMatrix22(int mp, double beta, double &dddcosb)
  {
    if((unsigned)abs(mp) > 2) throw std::domain_error("l >= |m'| is required.");
    if(mp == 2)
    {
      // d    = cos(b/2)^4
      // dddb = - 2 * sin(b/2) * cos(b/2)^3 = - sinb * cos(b/2)^2
      // dddcosb = - dddb / sinb = cos(b/2)^2
      auto cos_hb = cos(0.5 * beta);
      dddcosb = cos_hb * cos_hb;
      return dddcosb * dddcosb;
    }    
    else if(mp == 1)
    {
      // d    = sinb * cos(b/2)^2
      // dddb = - sinb^2 + cos(b/2)^2
      // dddcosb = - dddb / sinb = sinb - 1 / (2 * tan(b/2))
      auto sinb = sin(beta);
      auto cos_hb = cos(0.5 * beta);
      auto tan_hb = tan(0.5 * beta);
      dddcosb = sinb - 0.5 / tan_hb;
      return sinb * cos_hb * cos_hb;
    }    
    else if(mp == 0)
    {
      // sqrt(6) / 2 = 1.224744871391589
      // sqrt(6) / 4 = 0.6123724356957945
      // d    = sqrt(6) / 4 * sinb * sinb
      // dddb = sqrt(6) / 2 * sinb * cosb
      // dddcosb = - dddb / sinb = - sqrt(6) / 2 * cosb
      auto cosb = cos(beta);
      dddcosb = - 1.224744871391589 * cosb;
      return 0.6123724356957945 * (1. - cosb * cosb);
    }    
    else if(mp == -1)
    {
      // d    = sinb * sin(b/2)^2
      // dddb = sinb^2 - sin(b/2)^2
      // dddcosb = - dddb / sinb = - sinb + tan(b/2) / 2
      auto sinb = sin(beta);
      auto sin_hb = sin(0.5 * beta);
      auto tan_hb = tan(0.5 * beta);
      dddcosb = - sinb + 0.5 * tan_hb;
      return sinb * sin_hb * sin_hb;
    }    
    else // if(mp == -2)
    {
      // d    = sin(b/2)^4
      // dddb = 2 * cos(b/2) * sin(b/2)^3 = sinb * sin(b/2)^2
      // dddcosb = - dddb / sinb = - sin(b/2)^2
      auto sin_hb = sin(0.5 * beta);
      dddcosb = - sin_hb * sin_hb;
      return dddcosb * dddcosb;
    }
  }

  double discreteDelta(double w, double x)
  {
    const auto a = fabs(x / w);
    if(a < 1) return (1. - a) / w;
    else return 0;
  }

  double HeavisideStep(double x)
  {
    if(x > 0) return 1.;
    else return 0;
  }
}