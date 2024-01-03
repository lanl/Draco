//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   special_functions/F4.cc
 * \author Kent Budge
 * \date   Tue Sep 21 09:20:10 2004
 * \brief  Implementation of F4
 * \note   Copyright (C) 2010-2024 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "F4.hh"
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
 * This is a translation of an implementation from Chris Fryer, LANL T-6.  It is based on
 * approximations from Takahashi, El Eid, Hillebrandt, 1978 AA, 67, 185.  Like the Antlia
 * approximations, these are valid over the entire real axis, but do not match the 1e-12 accuracy
 * typical of Antlia.
 *
 * \param eta Dimensionless chemical potential \f$\eta\f$
 *
 * \return Value of \f$F_4(\eta)\f$
 *
 * \post \c Result>=0
 */
double F4(double const eta) {
  double f4;

  //--case where eta > 1e-3

  if (eta > 1.e-3) {
    double const eta2 = eta * eta;
    double const eta3 = eta * eta2;
    double const eta4 = eta * eta3;
    double const eta5 = eta * eta4;
    //        double const eta6=eta*eta5;
    if (eta <= 30.) {
      f4 = (0.2 * eta5 + 6.5797 * eta3 + 45.4576 * eta) / (1. - exp(-1.9484 * eta));
    } else {
      f4 = 0.2 * eta5 + 6.5797 * eta3 + 45.4576 * eta;
    }
  } else {
    double const expeta = exp(eta);
    f4 = 24. * expeta / (1. + 0.0287 * exp(0.9257 * eta));
  }
  return f4;
}

} // end namespace rtt_sf

//------------------------------------------------------------------------------------------------//
// end of F4.cc
//------------------------------------------------------------------------------------------------//
