//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   special_functions/ExpInt.hh
 * \author Paul Talbot
 * \date   Tue Jul 26 14:48:13 MDT 2011
 * \brief  Declare the ExpInt function templates.
 * \note   Copyright (C) 2011-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef special_functions_ExpInt_hh
#define special_functions_ExpInt_hh

#include "ds++/config.h"

namespace rtt_sf {

//! Compute general exponential integral, order n, argument x \f$ E_n(x) \f$.
double En(unsigned const n, double const x);

//! Compute exponential integral, argument x \f$ Ei(x) \f$.
double Ei(double const x);
} // namespace rtt_sf

#endif //special_functions_ExpInt

//------------------------------------------------------------------------------------------------//
// end of ExpInt.hh
//------------------------------------------------------------------------------------------------//
