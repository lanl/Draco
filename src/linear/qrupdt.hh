//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   linear/qrupdt.hh
 * \author Kent Budge
 * \brief  Update the QR decomposition of a square matrix
 * \note   Copyright (C) 2010-2026 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef linear_qrupdt_hh
#define linear_qrupdt_hh

namespace rtt_linear {
//! Update the QR decomposition of a square matrix.
template <class RandomContainer>
void qrupdt(RandomContainer &r, RandomContainer &qt, const unsigned n, RandomContainer &u,
            RandomContainer &v);

} // end namespace rtt_linear

#endif // linear_qrupdt_hh

//------------------------------------------------------------------------------------------------//
// end of linear/qrupdt.hh
//------------------------------------------------------------------------------------------------//
