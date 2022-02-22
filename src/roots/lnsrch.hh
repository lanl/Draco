//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   roots/lnsrch.hh
 * \author Kent Budge
 * \date   Tue Aug 10 13:21:58 2004
 * \brief  Reduce norm of a set of functions on a ray.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef roots_lnsrch_hh
#define roots_lnsrch_hh

namespace rtt_roots {

//! Reduce the norm of a set of functions on a ray.

template <class RandomContainer, class Function_N_to_N>
void lnsrch(const RandomContainer &xold, double fold, const RandomContainer &g, RandomContainer &p,
            RandomContainer &x, double &f, bool &check, RandomContainer &fvec,
            const Function_N_to_N &vecfunc, double ALF, double min_lambda);

} // end namespace rtt_roots

#include "lnsrch.i.hh"

#endif // roots_lnsrch_hh

//------------------------------------------------------------------------------------------------//
// end of roots/lnsrch.hh
//------------------------------------------------------------------------------------------------//
