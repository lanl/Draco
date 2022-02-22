//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   linear/svbksb.hh
 * \author Kent Budge
 * \date   Tue Aug 10 13:08:03 2004
 * \brief  Solve a linear system from its singular value decomposition.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef linear_svbksb_hh
#define linear_svbksb_hh

namespace rtt_linear {
//! Solve a linear system given its singular value decomposition.
template <class RandomContainer>
void svbksb(const RandomContainer &u, const RandomContainer &w, const RandomContainer &v,
            const unsigned m, const unsigned n, const RandomContainer &b, RandomContainer &x);

} // end namespace rtt_linear

#endif // linear_svbksb_hh

//------------------------------------------------------------------------------------------------//
// end of linear/svbksb.hh
//------------------------------------------------------------------------------------------------//
