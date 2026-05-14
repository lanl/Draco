//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   special_functions/F12.cc
 * \author Kent Budge
 * \brief  Implementation of F12
 * \note   Copyright (C) 2010-2026 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "F12.hh"
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
 * \return Value of \f$F_{1/2}(\eta)\f$
 *
 * \post \c Result>=0
 */
double F12(double x) {
  // coefficients of the expansion
  //  const double an = 0.5;
  constexpr int m1 = 7;
  constexpr int k1 = 7;
  constexpr int m2 = 10;
  constexpr int k2 = 11;

  const std::array<double, m1 + 1> a1 = {
      5.75834152995465e6, 1.30964880355883e7, 1.07608632249013e7, 3.93536421893014e6,
      6.42493233715640e5, 4.16031909245777e4, 7.77238678539648e2, 1.0};
  const std::array<double, k1 + 1> b1 = {6.49759261942269e6, 1.70750501625775e7, 1.69288134856160e7,
                                         7.95192647756086e6, 1.83167424554505e6, 1.95155948326832e5,
                                         8.17922106644547e3, 9.02129136642157e1};
  const std::array<double, m2 + 1> a2 = {4.85378381173415e-14,
                                         1.64429113030738e-11,
                                         3.76794942277806e-9,
                                         4.69233883900644e-7,
                                         3.40679845803144e-5,
                                         1.32212995937796e-3,
                                         2.60768398973913e-2,
                                         2.48653216266227e-1,
                                         1.08037861921488e0,
                                         1.91247528779676e0,
                                         1.0};
  const std::array<double, k2 + 1> b2 = {
      7.28067571760518e-14, 2.45745452167585e-11, 5.62152894375277e-9, 6.96888634549649e-7,
      5.02360015186394e-5,  1.92040136756592e-3,  3.66887808002874e-2, 3.24095226486468e-1,
      1.16434871200131e0,   1.34981244060549e0,   2.01311836975930e-1, -2.14562434782759e-2};

  if (x < 2.0) {
    double const xx = std::exp(x);
    double rn = xx + a1[m1 - 1];
    for (int i = m1 - 2; i >= 0; i--) {
      rn = rn * xx + a1[i];
    }
    double den = b1[k1];
    for (int i = k1 - 1; i >= 0; i--) {
      den = den * xx + b1[i];
    }
    return xx * rn / den;
  } else {
    double const xx = 1.0 / (x * x);
    double rn = xx + a2[m2 - 1];
    for (int i = m2 - 2; i >= 0; i--) {
      rn = rn * xx + a2[i];
    }
    double den = b2[k2];
    for (int i = k2 - 1; i >= 0; i--) {
      den = den * xx + b2[i];
    }
    return x * std::sqrt(x) * rn / den;
  }
}

} // end namespace rtt_sf

//------------------------------------------------------------------------------------------------//
// end of F12.cc
//------------------------------------------------------------------------------------------------//
