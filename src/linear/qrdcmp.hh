//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   linear/qrdcmp.hh
 * \author Kent Budge
 * \date   Mon Aug  9 13:17:31 2004
 * \brief  Calculate the Q-R decomposition of a square matrix.
 * \note   Copyright (C) 2016-2020 Triad National Security, LLC.
 *         All rights reserved.
 *
 * Long description.
 */
//------------------------------------------------------------------------------------------------//

//------------------------------------------------------------------------------------------------//

#ifndef linear_qrdcmp_hh
#define linear_qrdcmp_hh

namespace rtt_linear {
//! Compute the QR decomposition of a square matrix.
template <class RandomContainer>
bool qrdcmp(RandomContainer &a, unsigned n, RandomContainer &c,
            RandomContainer &d);

} // end namespace rtt_linear

#endif // linear_qrdcmp_hh

//------------------------------------------------------------------------------------------------//
// end of linear/qrdcmp.hh
//------------------------------------------------------------------------------------------------//
