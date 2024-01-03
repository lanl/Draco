//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   special_functions/F2inv.cc
 * \author Kent Budge
 * \date   Tue Sep 21 09:20:10 2004
 * \brief  Implementation of F2inv
 * \note   Copyright (C) 2010-2024 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "F2inv.hh"
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
 * Center.  This uses a rational function expansion to get the inverse of the fermi-dirac
 * integral. Reference: antia apjs 84,101 1993
 *
 * \param f Value of \f$F_2(\eta)\f$
 *
 * \pre \c f>0
 *
 * \return Dimensionless chemical potential \f$\eta\f$
 */
double F2inv(double f) {
  double constexpr an = 2.0;
  int constexpr m1 = 4;
  int constexpr k1 = 3;
  int constexpr m2 = 4;
  int constexpr k2 = 3;

  std::array<double, m1 + 1> constexpr a1 = {125.829, 35974.8, 7993.89, 307.849, 1.0e0};
  std::array<double, k1 + 1> constexpr b1 = {251.657, 71940.5, 11494.2, -0.0140884};
  std::array<double, m2 + 1> constexpr a2 = {-3.04879, 714344, 23834.7, -1.08562e6, 1.0e0};
  std::array<double, k2 + 1> constexpr b2 = {-9.14637, 495569, 13194.3, 48419.5};

  if (f < 4.0e0) {
    double rn = f + a1[m1 - 1];
    for (int i = m1 - 2; i >= 0; i--) {
      rn = rn * f + a1[i];
    }
    double den = b1[k1];
    for (int i = k1 - 1; i >= 0; i--) {
      den = den * f + b1[i];
    }
    return log(f * rn / den);
  } else {
    double ff = 1.0 / std::pow(f, (1.0 / (1.0 + an)));
    double rn = ff + a2[m2 - 1];
    for (int i = m2 - 2; i >= 0; i--) {
      rn = rn * ff + a2[i];
    }
    double den = b2[k2];
    for (int i = k2 - 1; i >= 0; i--) {
      den = den * ff + b2[i];
    }
    return rn / (den * ff);
  }
}

//------------------------------------------------------------------------------------------------//
/*!
 * The Fermi-Dirac integral is defined as
 * \f[
 * F_n(\eta) = \int_0^\infty \frac{x^n}{e^{x-\eta}+1} dx
 * \f]
 *
 * This implementation is a translation of an implementation from the Chicago Astrophysical Flash
 * Center.  This uses a rational function expansion to get the inverse of the fermi-dirac
 * integral. Reference: antia apjs 84,101 1993
 *
 * \param[in] f Value of \f$F_2(\eta)\f$
 * \param[out] eta Dimensionless chemical potential \f$\eta\f$
 * \param[out] deta Derivative of dimensionless chemical potential
 * \f$\frac{d\eta}{dF_2(\eta)}\f$
 *
 * \pre \c f>0
 */
void F2inv(double const f, double &eta, double &deta) {
  double constexpr an = 2.0;
  int constexpr m1 = 4;
  int constexpr k1 = 3;
  int constexpr m2 = 4;
  int constexpr k2 = 3;
  std::array<double, m1 + 1> constexpr a1 = {125.829, 35974.8, 7993.89, 307.849, 1.0e0};
  std::array<double, k1 + 1> constexpr b1 = {251.657, 71940.5, 11494.2, -0.0140884};
  std::array<double, m2 + 1> constexpr a2 = {-3.04879, 714344, 23834.7, -1.08562e6, 1.0e0};
  std::array<double, k2 + 1> constexpr b2 = {-9.14637, 495569, 13194.3, 48419.5};

  if (f < 4.0e0) {
    double rn = f + a1[m1 - 1];
    double drndf = 1;
    for (int i = m1 - 2; i >= 0; i--) {
      drndf = drndf * f + rn;
      rn = rn * f + a1[i];
    }
    double den = b1[k1];
    double ddendf = 0;
    for (int i = k1 - 1; i >= 0; i--) {
      ddendf = ddendf * f + den;
      den = den * f + b1[i];
    }
    eta = log(f * rn / den);
    deta = (rn + f * drndf - f * rn * ddendf / den) / (f * rn);
  } else {
    double ff = 1.0 / std::pow(f, (1.0 / (1.0 + an)));
    double dffdf = -(1.0 / (1.0 + an)) * ff / f;
    double rn = ff + a2[m2 - 1];
    double drndf = dffdf;
    for (int i = m2 - 2; i >= 0; i--) {
      drndf = drndf * ff + rn * dffdf;
      rn = rn * ff + a2[i];
    }
    double den = b2[k2];
    double ddendf = 0;
    for (int i = k2 - 1; i >= 0; i--) {
      ddendf = ddendf * ff + den * dffdf;
      den = den * ff + b2[i];
    }
    eta = rn / (den * ff);
    deta = drndf / (den * ff) - rn * (ddendf * ff + den * dffdf) / (den * den * ff * ff);
  }
}

} // end namespace rtt_sf

//------------------------------------------------------------------------------------------------//
// end of F2inv.cc
//------------------------------------------------------------------------------------------------//
