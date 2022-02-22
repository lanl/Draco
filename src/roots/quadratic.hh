//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   roots/quadratic.hh
 * \author Kent Budge
 * \date   Wed Aug 11 15:21:38 2004
 * \brief  Solve a quadratic equation.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef roots_quadratic_hh
#define roots_quadratic_hh

#include "ds++/Assert.hh"
#include <cmath>

namespace rtt_roots {

//------------------------------------------------------------------------------------------------//
/*! Solve a quadratic equation.
 *
 * \arg \a Field A real field type, such as float or double.
 *
 * \param[in] a
 * Quadratic coefficient
 * \param[in] b
 * Linear[in] coefficient
 * \param[in] c
 * Constant coefficient
 * \param[out] r1
 * First root
 * \param[out] r2
 * Second root
 *
 * \note The roots are not returned in any particular order.  The client must decide which root to
 * use for a particular application.
 */

template <class Field>
void quadratic(Field const &a, Field const &b, Field const &c, Field &r1, Field &r2) {
  using namespace std;

  Field det = sqrt(b * b - 4. * a * c);
  if (b < 0.0)
    det = -det;
  Field const q = -0.5 * (b + det);
  r1 = q / a;
  r2 = c / q;
}

} // end namespace rtt_roots

#endif // roots_quadratic_hh

//------------------------------------------------------------------------------------------------//
// end of roots/quadratic.hh
//------------------------------------------------------------------------------------------------//
