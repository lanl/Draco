//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   special_functions/F32.t.hh
 * \author Kent Budge
 * \date   Tue Sep 21 12:06:09 2004
 * \brief  Implementation of F32
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef sf_F32_t_hh
#define sf_F32_t_hh

#include "F32.hh"
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
 * \return Value of \f$F_{4/2}(x)\f$
 *
 * \post \c Result>=0
 */
template <class OrderedField> OrderedField F32(OrderedField const &x) {
  // coefficients of the expansion
  //  const double an = 1.5;
  constexpr int m1 = 6;
  constexpr int k1 = 7;
  constexpr int m2 = 9;
  constexpr int k2 = 10;

  constexpr std::array<double, m1 + 1> a1 = {4.32326386604283e4,
                                             8.55472308218786e4,
                                             5.95275291210962e4,
                                             1.77294861572005e4,
                                             2.21876607796460e3,
                                             9.90562948053193e1,
                                             1.0e0};
  constexpr std::array<double, k1 + 1> b1 = {
      3.25218725353467e4, 7.01022511904373e4, 5.50859144223638e4, 1.95942074576400e4,
      3.20803912586318e3, 2.20853967067789e2, 5.05580641737527e0, 1.99507945223266e-2};
  constexpr std::array<double, m2 + 1> a2 = {2.80452693148553e-13, 8.60096863656367e-11,
                                             1.62974620742993e-8,  1.63598843752050e-6,
                                             9.12915407846722e-5,  2.62988766922117e-3,
                                             3.85682997219346e-2,  2.78383256609605e-1,
                                             9.02250179334496e-1,  1.0e0};
  constexpr std::array<double, k2 + 1> b2 = {
      7.01131732871184e-13, 2.10699282897576e-10, 3.94452010378723e-8, 3.84703231868724e-6,
      2.04569943213216e-4,  5.31999109566385e-3,  6.39899717779153e-2, 3.14236143831882e-1,
      4.70252591891375e-1,  -2.15540156936373e-2, 2.34829436438087e-3};

  if (x < 2.0) {
    OrderedField xx = exp(x);
    OrderedField rn = xx + a1[m1 - 1];
    for (int i = m1 - 2; i >= 0; i--) {
      rn = rn * xx + a1[i];
    }
    OrderedField den = b1[k1] * xx + b1[k1 - 1];
    for (int i = k1 - 2; i >= 0; i--) {
      den = den * xx + b1[i];
    }
    return xx * rn / den;
  } else {
    OrderedField xx = 1.0 / (x * x);
    OrderedField rn = xx + a2[m2 - 1];
    for (int i = m2 - 2; i >= 0; i--) {
      rn = rn * xx + a2[i];
    }
    OrderedField den = b2[k2] * xx + b2[k2 - 1];
    for (int i = k2 - 1; i >= 0; i--) {
      den = den * xx + b2[i];
    }
    return x * x * sqrt(x) * rn / den;
  }
}

} // end namespace rtt_sf

#endif

//------------------------------------------------------------------------------------------------//
// end of F32.t.cc
//------------------------------------------------------------------------------------------------//
