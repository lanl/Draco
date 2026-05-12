//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   roots/quadratic.cc
 * \author Kent Budge
 * \brief  Specializations of quadratic
 * \note   Copyright (C) 2010-2026 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "quadratic.hh" // NOLINT(misc-include-cleaner)

namespace rtt_roots {

template void quadratic(double const &a, double const &b, double const &c, double &r1, double &r2);

} // end namespace rtt_roots

//------------------------------------------------------------------------------------------------//
// end of quadratic.cc
//------------------------------------------------------------------------------------------------//
