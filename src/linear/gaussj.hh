//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   linear/gaussj.hh
 * \author Kent Budge
 * \brief  Solve a system of equations by Gaussian elimination.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef linear_gaussj_hh
#define linear_gaussj_hh

namespace rtt_linear {

//------------------------------------------------------------------------------------------------//
//! Solve a system of linear equations.  Single-subscript computed-index version.
template <class RandomContainer>
void gaussj(RandomContainer &A, unsigned n, RandomContainer &b, unsigned m);

//! Solve a system of linear equations. Double-subscript version.
template <class DoubleRandomContainer, class RandomContainer>
void gaussj(DoubleRandomContainer &A, RandomContainer &b);

} // end namespace rtt_linear

#endif // linear_gaussj_hh

//------------------------------------------------------------------------------------------------//
// end of linear/gaussj.hh
//------------------------------------------------------------------------------------------------//
