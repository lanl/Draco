//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   norms/norm.hh
 * \author Kent Budge
 * \date   Tue Sep 18 08:53:41 2007
 * \brief  Declare template class norm
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef norms_norm_hh
#define norms_norm_hh

#include "ds++/DracoMath.hh"

namespace rtt_norms {

//------------------------------------------------------------------------------------------------//
//! Compute the norm of a scalar.
template <typename Field> double norm(Field const &x) { return x * x; }

//------------------------------------------------------------------------------------------------//
//! Compute the norm of the difference of two scalars.
template <typename Field1, typename Field2> double norm_diff(Field1 const &x, Field2 const &y) {
  return rtt_dsxx::square(x - y);
}

} // end namespace rtt_norms

#endif // norms_norm_hh

//------------------------------------------------------------------------------------------------//
// end of norms/norm.hh
//------------------------------------------------------------------------------------------------//
