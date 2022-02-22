//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   linear/svdcmp.hh
 * \author Kent Budge
 * \date   Mon Aug  9 13:17:31 2004
 * \brief  Calculate the singular value decomposition of a matrix.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef linear_svdcmp_hh
#define linear_svdcmp_hh

namespace rtt_linear {

//! Compute the singular value decomposition of a matrix.
template <class RandomContainer>
void svdcmp(RandomContainer &a, const unsigned m, const unsigned n, RandomContainer &w,
            RandomContainer &v);

} // end namespace rtt_linear

#endif // linear_svdcmp_hh

//------------------------------------------------------------------------------------------------//
// end of linear/svdcmp.hh
//------------------------------------------------------------------------------------------------//
