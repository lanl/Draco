//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   special_functions/F_eta.cc
 * \author Kent Budge
 * \date   Mon Sep 20 15:01:53 2004
 * \brief  Implementation of F_eta. These routines are based on C routines from Numerical Recipes.
 * \note   Copyright (C) 2010-2024 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "F_eta.hh"
#include "Factorial.hh"
#include "ds++/DracoMath.hh"
#include "ode/quad.hh"
#include "ode/rkqs.hh"
#include "units/PhysicalConstants.hh"

namespace rtt_sf {
using namespace std;
using namespace rtt_dsxx;
using namespace rtt_ode;
using rtt_units::PI;

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Evaluate the relativistic Fermi-Dirac integral
 *
 * The relativistic Fermi-Dirac integral is defined as
 * \f[
 * F_{3/2}(\eta, \gamma) = \frac{1}{2^{3/2}\gamma^{5/2}} \int_0^\infty
 * \frac{(x^2+2x)^{3/2}}{e^\frac{x-\eta}{\gamma}+1} dx
 * \f]
 * The dimensionless number density is its partial derivative with eta.
 *
 * \param[in] eta Dimensionless chemical potential \f$\eta=\frac{\mu}{kT}\f$
 * \param[in] gamma Dimensionless temperature \f$\gamma=\frac{kT}{mc^2}\f$
 *
 * \return Dimensionless number density \f$\frac{\partial F_{3/2}}{\partial\eta}=\frac{3Nh^3}{8\pi
 *         m^3c^3}\f$.
 */
double F_eta(double const eta, double const gamma) {
  Require(gamma > 0.0);

  double const TOL = 1.0e-5;

  double const e = exp(eta / gamma);
  double const de = e / gamma;
  if (e <= 0.5 && gamma < 0.2) {
    // Classical regime.  Expand in powers of exp(eta/gamma).
    double sum = 0.0;
    double dsum = 0.0;
    double fac = 1;
    for (int i = 0; i < 9; i++) {
      double si = 1;
      double dsi = 0;
      double ep = 1;
      double dep = 0;
      double sign = 1;
      for (int j = 2;; j++) {
        dep = dep * e + ep * dep;
        ep *= e;
        sign *= -1;
        double srt = pow(static_cast<double>(j), i + 1.5);
        double term = sign * ep / srt;
        if (fabs(term) < fabs(si) * std::numeric_limits<double>::epsilon())
          break;
        double dterm = sign * dep / srt;
        si += term;
        dsi += dterm;
      }
      double const term = 0.75 * sqrt(PI) * e * fac * si / factorial(i);
      if (fabs(term) < fabs(sum) * std::numeric_limits<double>::epsilon())
        break;
      double const dterm = 0.75 * sqrt(PI) * fac * (de * si + e * dsi) / factorial(i);
      sum += term;
      dsum += dterm;
      fac *= 0.25 * gamma * (3 - 2 * i) * (2.5 + i);
    }
    return dsum;
  } else {
    // Degenerate regime?
    double const ee = square(eta + 1) - 1;

    if (ee <= 0) {
      Insist(false, std::string("Please add a unit test for this case and ") +
                        "then re-enable Feta_brute(eta,gamma).");
      // return Feta_brute(eta, gamma);
    } else {
      double const dde = 2 * (eta + 1);
      double const x = sqrt(ee);
      double const dx = 0.5 * dde / x;
      double const rad = sqrt(x * x + 1);
      double const drad = x * dx / rad;

      double n1, dn1;
      // Two versions of Sommerfeld expansion with different roundoff properties.
      if (x < 0.07) {
        n1 = 0.2 * x * x * x * x * x *
             (1 - (5. / 14.) * x * x + (5. / 24.) * x * x * x * x -
              (25. / 176.) * x * x * x * x * x * x + (35. / 1024.) * x * x * x * x * x * x * x * x);
        dn1 = dx * (5 * n1 / x + 0.2 * x * x * x * x * x * x *
                                     (-5. / 7. + (5. / 6.) * x * x - (75. / 88.) * x * x * x * x +
                                      (35. / 128.) * x * x * x * x * x * x));
      } else {
        n1 = (x * (2 * x * x - 3) * rad + 3 * log(rad + x)) / 8;
        dn1 = (dx * (2 * x * x - 3) * rad + 4 * x * x * dx * rad + x * (2 * x * x - 3) * drad +
               3 * (drad + dx) / (rad + x)) /
              8;
      }
      double const n3 = 7 * square(square(PI * gamma)) * (2 * x * x - 1) * rad / (120 * x * x * x);
      double const dn3 = 7 * square(square(PI * gamma / x)) * dx / (40 * rad);
      if (fabs(n3 / n1) < TOL) {
        // Sommerfeld expansion is sufficiently accurate
        double const dn2 = square(PI * gamma) * (2 * x * x + 1) * dx / (2 * rad);
        return (dn1 + dn2 + dn3) / (2 * gamma * gamma * sqrt(2 * gamma));
      } else {
        Insist(false, std::string("Please add a unit test for this case and ") +
                          "then re-enable Feta_brute(eta,gamma).");
        // return Feta_brute(eta, gamma);
      }
    }
  }
}

} // end namespace rtt_sf

//------------------------------------------------------------------------------------------------//
// end of F_eta.cc
//------------------------------------------------------------------------------------------------//
