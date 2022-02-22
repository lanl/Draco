//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   roots/cubic1.hh
 * \author Kent Budge
 * \date   Wed Sep 15 10:04:02 MDT 2010
 * \brief  Solve a cubic equation assumed to have one real root
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef roots_cubic1_hh
#define roots_cubic1_hh

#include "ds++/config.h"

namespace rtt_roots {

//------------------------------------------------------------------------------------------------//
/*! Solver for cubic equation known to have only one real root.
 *
 *
 * \arg \a Field A class representing the real numbers
 *
 * The parameters correspond to the coefficients in a cubic equation \f$ x^3 + ax^2 + bx + c = 0\f$
 */
template <typename Field> Field cubic1(Field const &a, Field const &b, Field const &c);

} // end namespace rtt_roots

#endif
// roots_cubic1_hh

//------------------------------------------------------------------------------------------------//
// end of cubic1.hh
//------------------------------------------------------------------------------------------------//
