//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   special_functions/F2.cc
 * \author Kent Budge
 * \brief  Implementation of F2
 * \note   Copyright (C) 2010-2026 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "F2.hh"
#include <array>
#include <cmath>

namespace rtt_sf {
using namespace std;

//------------------------------------------------------------------------------------------------//
/*!
 * The Fermi-Dirac integral is defined as
 * \f[
 * F_n(\eta) = \int_0^\infty \frac{x^n}{e^{x-\eta}+1} dx
 * \f]
 *
 * This implementation is a translation of an implementation from the Chicago Astrophysical Flash
 * Center.  This uses a rational function expansion to get the fermi-dirac integral. Reference:
 * antia apjs 84,101 1993
 *
 * \param x Dimensionless chemical potential \f$\eta\f$
 *
 * \return Value of \f$F_2(\eta)\f$
 *
 * \post \c Result>=0
 */
double F2(double const x) {
  //..load the coefficients of the expansion
  unsigned constexpr m1 = 7;
  unsigned constexpr k1 = 4;
  unsigned constexpr m2 = 5;
  unsigned constexpr k2 = 9;

  std::array<double, m1 + 1> constexpr a1 = {
      -1.434885992395e8, -2.001711155617e8, -8.507067153428e7, -1.175118281976e7,
      -3.145120854293e5, 4.275771034579e3,  -8.069902926891e1, 1.0e0};
  std::array<double, k1 + 1> constexpr b1 = {-7.174429962316e7, -1.090535948744e8,
                                             -5.350984486022e7, -9.646265123816e6,
                                             -5.113415562845e5};
  std::array<double, m2 + 1> constexpr a2 = {6.919705180051e-8, 1.134026972699e-5,
                                             7.967092675369e-4, 2.432500578301e-2,
                                             2.784751844942e-1, 1.0e0};
  std::array<double, k2 + 1> constexpr b2 = {
      2.075911553728e-7,  3.197196691324e-5,  2.074576609543e-3, 5.250009686722e-2,
      3.171705130118e-1,  -1.147237720706e-1, 6.638430718056e-2, -1.356814647640e-2,
      -3.648576227388e-2, 3.621098757460e-2};

  if (x < 2.0e0) {
    double const xx = exp(x);
    double rn = xx + a1[m1 - 1];
    for (int i = m1 - 2; i >= 0; --i) {
      rn = rn * xx + a1[i];
    }
    double den = b1[k1];
    for (int i = k1 - 1; i >= 0; --i) {
      den = den * xx + b1[i];
    }
    double const Result = xx * rn / den;
    return Result;
  } else {
    double const xx = 1.0 / (x * x);
    double rn = xx + a2[m2 - 1];
    for (int i = m2 - 2; i >= 0; --i) {
      rn = rn * xx + a2[i];
    }
    double den = b2[k2];
    for (int i = k2 - 1; i >= 0; --i) {
      den = den * xx + b2[i];
    }
    double const Result = x * x * x * rn / den;
    return Result;
  }
}

} // end namespace rtt_sf

//------------------------------------------------------------------------------------------------//
// end of F2.cc
//------------------------------------------------------------------------------------------------//
