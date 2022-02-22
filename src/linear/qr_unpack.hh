//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   linear/qr_unpack.hh
 * \author Kent Budge
 * \date   Mon Aug  9 13:17:31 2004
 * \brief  Compute an explicit representation of a packed QR decomposition.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef linear_qr_unpack_hh
#define linear_qr_unpack_hh

namespace rtt_linear {
//! Compute an explicit representation of a packed QR decomposition.
template <class RandomContainer>
void qr_unpack(RandomContainer &r, const unsigned n, const RandomContainer &c,
               const RandomContainer &d, RandomContainer &qt);

} // end namespace rtt_linear

#endif // linear_qr_unpack_hh

//------------------------------------------------------------------------------------------------//
// end of linear/qr_unpack.hh
//------------------------------------------------------------------------------------------------//
