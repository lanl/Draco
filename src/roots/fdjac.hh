//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   roots/fdjac.hh
 * \author Kent Budge
 * \date   Wed Aug 11 08:07:04 2004
 * \brief  Compute the Jacobian of a nonlinear system of equations
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef roots_fdjac_hh
#define roots_fdjac_hh

#include "ds++/Assert.hh"
#include "ds++/Soft_Equivalence.hh"
#include <algorithm>
#include <cfloat>
#include <cmath>
#include <limits>
#include <vector>

namespace rtt_roots {

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Calculate the Jacobian of a nonlinear system of equations.
 *
 * This procedure computes the Jacobian using a forward-difference approximation.
 *
 * \arg \a Field A field type
 * \arg \a Function_N_to_N A function representing a set of N functions of N variables.
 *
 * \param x Point at which the Jacobian is to be evaluated.
 * \param fvec Residuals of the equations at x.
 * \param df On return, contains the Jacobian. The ordering is that df[i+n*j] contains the jth
 *         derivative of the ith residual.
 * \param vecfunc Multifunctor returning the residuals of the nonlinear equations.
 *
 * \pre \c x.size()==fvec.size()
 * \post \c df.size()==square(x.size())
 */

template <class Field, class Function_N_to_N>
void fdjac(const std::vector<Field> &x, const std::vector<Field> &fvec, std::vector<Field> &df,
           const Function_N_to_N &vecfunc) {
  Require(x.size() == fvec.size());

  using std::abs;
  using std::numeric_limits;
  using std::vector;

  // Square root of the machine precision
  static const double EPS = sqrt(numeric_limits<Field>::epsilon());

  const size_t n = x.size();

  df.resize(n * n);

  vector<Field> f(n);
  vector<Field> xt = x;
  for (unsigned j = 0; j < n; j++) {
    Field temp = xt[j];
    Field h = EPS * abs(temp);
    if (std::abs(h) < static_cast<Field>(FLT_MIN))
      h = EPS;
    xt[j] = temp + h;
    h = xt[j] - temp;
    vecfunc(xt, f);
    xt[j] = temp;
    for (unsigned i = 0; i < n; i++) {
      df[i + n * j] = (f[i] - fvec[i]) / h;
    }
  }
}

} // end namespace rtt_roots

#endif // roots_fdjac_hh

//------------------------------------------------------------------------------------------------//
// end of roots/fdjac.hh
//------------------------------------------------------------------------------------------------//
