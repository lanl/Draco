//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   roots/zbrent.hh
 * \author Kent Budge
 * \date   Tue Aug 17 15:57:06 2004
 * \brief  Find a bracketed root of a function.
 * \note   Copyright (C) 2010-2024 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef roots_zbrent_hh
#define roots_zbrent_hh

namespace rtt_roots {

//------------------------------------------------------------------------------------------------//
// Use explicit instantiation for these functions.
//------------------------------------------------------------------------------------------------//

//! Pinpoint a bracketed root of a function.
template <typename Function, typename Real>
Real zbrent(Function func, Real x1, Real x2, unsigned itmax, Real &tol, Real &ftol);

} // end namespace rtt_roots

// Use implicit template instantiation.
#include "zbrent.i.hh"

#endif // roots_zbrent_hh

//------------------------------------------------------------------------------------------------//
// end of roots/zbrent.hh
//------------------------------------------------------------------------------------------------//
