//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   norms/L2norm.hh
 * \author Kent Budge
 * \date   Tue Sep 18 08:22:08 2007
 * \brief  Define template function L2norm
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef norms_L2norm_hh
#define norms_L2norm_hh

namespace rtt_norms {
//! Compute the L2-norm of a vector.
template <typename In> double L2norm(In const &x);

//! Compute the L2-norm of the difference between two vectors.
template <typename In1, typename In2> double L2norm_diff(In1 const &x, In2 const &y);
} // end namespace rtt_norms

#include <norms/L2norm.i.hh>

#endif // norms_L2norm_hh

//------------------------------------------------------------------------------------------------//
// end of norms/L2norm.hh
//------------------------------------------------------------------------------------------------//
