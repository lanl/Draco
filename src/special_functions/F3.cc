//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   special_functions/F3.cc
 * \author Kent Budge
 * \date   Tue Sep 21 09:20:10 2004
 * \brief  Implementation of F3
 * \note   Copyright (C) 2010-2024 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "F3.hh"
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
 * \return Value of \f$F_3(\eta)\f$
 *
 * \post \c Result>=0
 */
double F3(double const x) {

  //..load the coefficients of the expansion
  unsigned constexpr m1 = 4;
  unsigned constexpr k1 = 6;
  unsigned constexpr m2 = 7;
  unsigned constexpr k2 = 7;

  std::array<double, m1 + 1> constexpr a1 = {6.317036716422e2, 7.514163924637e2, 2.711961035750e2,
                                             3.274540902317e1, 1.0e0};

  std::array<double, k1 + 1> constexpr b1 = {
      1.052839452797e2,  1.318163114785e2,  5.213807524405e1,  7.500064111991e0,
      3.383020205492e-1, 2.342176749453e-3, -8.445226098359e-6};

  std::array<double, m2 + 1> constexpr a2 = {
      1.360999428425e-8, 1.651419468084e-6, 1.021455604288e-4, 3.041270709839e-3,
      4.584298418374e-2, 3.440523212512e-1, 1.077505444383e0,  1.0e0};

  std::array<double, k2 + 1> constexpr b2 = {
      5.443997714076e-8, 5.531075760054e-6, 2.969285281294e-4, 6.052488134435e-3,
      5.041144894964e-2, 1.048282487684e-1, 1.280969214096e-2, -2.851555446444e-3};

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
    double const Result = x * x * x * x * rn / den;
    return Result;
  }
}

} // end namespace rtt_sf

//------------------------------------------------------------------------------------------------//
// end of F3.cc
//------------------------------------------------------------------------------------------------//
