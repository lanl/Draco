//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   special_functions/F1.cc
 * \author Kent Budge
 * \brief  Implementation of F1
 * \note   Copyright (C) 2010-2024 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "F1.hh"
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
 * This implementation is modelled after Fermi-Dirac implementations from the Chicago Astrophysical
 * Flash Center.  These use a rational function expansion to get the fermi-dirac
 * integral. Reference: antia apjs 84,101 1993
 *
 * \param x Dimensionless chemical potential \f$\eta\f$
 * \return Value of \f$F_1(\eta)\f$
 *
 * \post \c Result>=0
 */
double F1(double const x) {

  //..load the coefficients of the expansion
  unsigned constexpr m1 = 7;
  unsigned constexpr k1 = 4;
  unsigned constexpr m2 = 9;
  unsigned constexpr k2 = 5;

  std::array<double, m1 + 1> const a1 = {
      -7.606458638543e7, -1.143519707857e8, -5.167289383236e7, -7.304766495775e6,
      -1.630563622280e5, 3.145920924780e3,  -7.156354090495e1, 1.0};

  std::array<double, k1 + 1> const b1 = {-7.606458639561e7, -1.333681162517e8, -7.656332234147e7,
                                         -1.638081306504e7, -1.044683266663e6};

  std::array<double, m2 + 1> const a2 = {-3.493105157219e-7, -5.628286279892e-5,
                                         -5.188757767899e-3, -2.097205947730e-1,
                                         -3.353243201574e0,  -1.682094530855e1,
                                         -2.042542575231e1,  3.551366939795e0,
                                         -2.400826804233e0,  1.0};

  std::array<double, k2 + 1> const b2 = {-6.986210315105e-7, -1.102673536040e-4, -1.001475250797e-2,
                                         -3.864923270059e-1, -5.435619477378e0,  -1.563274262745e1};

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
    double const Result = x * x * rn / den;
    return Result;
  }
}

} // end namespace rtt_sf

//------------------------------------------------------------------------------------------------//
// end of F1.cc
//------------------------------------------------------------------------------------------------//
