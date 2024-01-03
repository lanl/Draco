//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   special_functions/F_eta_inv.cc
 * \author Kent Budge
 * \date   Mon Sep 20 15:01:53 2004
 * \brief  Implementation of F_eta_inv.
 * \note   Copyright (C) 2010-2024 Triad National Security, LLC., All rights reserved.
 *
 * The implementations here are derived from C implementations from Numerical Recipes.*/
//------------------------------------------------------------------------------------------------//

#include "F_eta_inv.hh"
#include "F12inv.hh"
#include "F_eta.hh"
#include "ode/quad.hh"
#include "ode/rkqs.hh"
#include "roots/zbrac.hh"
#include "roots/zbrent.hh"

namespace rtt_sf {
using namespace std;
using namespace rtt_roots;

// Parametrization of integrand and inverse functions
static double lgamma, ln;

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Return residual of inversion of F_eta.
 *
 * \param eta Dimensionless chemical potential
  * \return Value of the relativistic Fermi-Dirac integral for the given chemical potential.
 */
static double Feta_diff(double eta) { return F_eta(eta, lgamma) - ln; }

//------------------------------------------------------------------------------------------------//
/*!
 * The relativistic Fermi-Dirac integral is defined as
 * \f[
 * F_{3/2}(\eta, \gamma) = \frac{1}{2^{3/2}\gamma^{5/2}} \int_0^\infty
 * \frac{(x^2+2x)^{3/2}}{e^\frac{x-\eta}{\gamma}+1} dx
 * \f]
 * The dimensionless number density is its partial derivative with eta.
 *
 * \param n Dimensionless number density \f$\frac{\partial F_{3/2}}{\partial\eta}=\frac{3Nh^3}{\pi g
 *       m^3c^3}\f$
 * \param gamma Dimensionless temperature \f$\gamma=\frac{kT}{mc^2}\f$
 *
 * \pre \c n>0
 * \pre \c gamma>0
 *
 * \return Dimensionless chemical potential \f$\eta=\frac{\mu}{kT}\f$
 *
 * \note This implementation is very expensive.  it is also of limited utility, since for problems
 *       for which the relativistic form of the Fermi-Dirac function is needed, pair production is
 *       likely to be important, requiring inversion of \f$F_{3/2}(\eta, \gamma) - F_{3/2}(-2-\eta,
 *       \gamma)\f$ rather than just \f$F_{3/2}(\eta, \gamma)\f$
 */
double F_eta_inv(double const n, double const gamma) {
  Require(n >= 0);
  Require(gamma > 0);

  if (std::abs(n) < numeric_limits<double>::min())
    return -numeric_limits<double>::max();

  const double TOL = 1.0e-5;

  lgamma = gamma;
  ln = n;

  // Estimate eta from cold limit
  double eta = gamma * F12inv(2 * n * gamma / 3);

  // Estimate eta from cold degenerate limit
  double etad = sqrt(1. + pow(n, 2. / 3.));
  eta = min(eta, etad);

  double x1 = eta * 0.99;
  double x2 = 1.01 * eta;
  if (x2 < x1)
    swap(x1, x2);
  zbrac(&Feta_diff, x1, x2);
  double tol = TOL * gamma;
  double ftol = 0;
  eta = zbrent(Feta_diff, x1, x2, 100, tol, ftol);

  return eta;
}

} // end namespace rtt_sf

//------------------------------------------------------------------------------------------------//
// end of F_eta_inv.cc
//------------------------------------------------------------------------------------------------//
