//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   special_functions/FM12.t.hh
 * \author Kent Budge
 * \date   Tue Sep 21 12:06:09 2004
 * \brief  Implementation of FM12
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef sf_FM12_t_hh
#define sf_FM12_t_hh

#include "FM12.hh"
#include <array>
#include <cmath>

namespace rtt_sf {

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
 * \return Value of \f$F_{4/2}(\eta)\f$
 *
 * \post \c Result>=0
 */
template <class OrderedField> OrderedField FM12(OrderedField const &x) {
  // coefficients of the expansion
  //  const double an = -0.5;
  constexpr int m1 = 7;
  constexpr int k1 = 7;
  constexpr int m2 = 11;
  constexpr int k2 = 11;

  constexpr std::array<double, m1 + 1> a1 = {
      1.71446374704454e7, 3.88148302324068e7, 3.16743385304962e7, 1.14587609192151e7,
      1.83696370756153e6, 1.14980998186874e5, 1.98276889924768e3, 1.0e0};
  constexpr std::array<double, k1 + 1> b1 = {
      9.67282587452899e6, 2.87386436731785e7, 3.26070130734158e7, 1.77657027846367e7,
      4.81648022267831e6, 6.13709569333207e5, 3.13595854332114e4, 4.35061725080755e2};
  constexpr std::array<double, m2 + 1> a2 = {
      -4.46620341924942e-15, -1.58654991146236e-12, -4.44467627042232e-10, -6.84738791621745e-8,
      -6.64932238528105e-6,  -3.69976170193942e-4,  -1.12295393687006e-2,  -1.60926102124442e-1,
      -8.52408612877447e-1,  -7.45519953763928e-1,  2.98435207466372e0,    1.0e0};
  constexpr std::array<double, k2 + 1> b2 = {
      -2.23310170962369e-15, -7.94193282071464e-13, -2.22564376956228e-10, -3.43299431079845e-8,
      -3.33919612678907e-6,  -1.86432212187088e-4,  -5.69764436880529e-3,  -8.34904593067194e-2,
      -4.78770844009440e-1,  -4.99759250374148e-1,  1.86795964993052e0,    4.16485970495288e-1};

  if (x < 2.0) {
    OrderedField xx = exp(x);
    OrderedField rn = xx + a1[m1 - 1];
    for (int i = m1 - 2; i >= 0; i--) {
      rn = rn * xx + a1[i];
    }
    OrderedField den = b1[k1];
    for (int i = k1 - 1; i >= 0; i--) {
      den = den * xx + b1[i];
    }
    return xx * rn / den;
  } else {
    OrderedField xx = 1.0 / (x * x);
    OrderedField rn = xx + a2[m2 - 1];
    for (int i = m2 - 2; i >= 0; i--) {
      rn = rn * xx + a2[i];
    }
    OrderedField den = b2[k2];
    for (int i = k2 - 1; i >= 0; i--) {
      den = den * xx + b2[i];
    }
    return sqrt(x) * rn / den;
  }
}

} // end namespace rtt_sf

#endif

//------------------------------------------------------------------------------------------------//
// end of FM12.t.cc
//------------------------------------------------------------------------------------------------//
