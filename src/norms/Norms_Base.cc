//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   Norms_Base.cc
 * \author Rob Lowrie
 * \date   Fri Jan 14 13:01:37 2005
 * \brief  Implemention for Norms_Base class.
 * \note   Copyright (C) 2005-2022 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

#include "Norms_Base.hh"
#include "ds++/Soft_Equivalence.hh"

using namespace rtt_norms;

//------------------------------------------------------------------------------------------------//
//! Default constructor.
Norms_Base::Norms_Base() : d_sum_L1(-42.0), d_sum_L2(-42.0), d_Linf(-42.0), d_sum_weights(-42.0) {
  Norms_Base::reset();
}

//------------------------------------------------------------------------------------------------//
//!  Re-initializes the norm values.
void Norms_Base::reset() {
  d_sum_L1 = 0.0;
  d_sum_L2 = 0.0;
  d_Linf = 0.0;
  d_sum_weights = 0.0;
}
//------------------------------------------------------------------------------------------------//
//! Equality operator.
bool Norms_Base::operator==(const Norms_Base &n) const {
  return rtt_dsxx::soft_equiv(d_sum_L1, n.d_sum_L1) && rtt_dsxx::soft_equiv(d_sum_L2, n.d_sum_L2) &&
         rtt_dsxx::soft_equiv(d_Linf, n.d_Linf) &&
         rtt_dsxx::soft_equiv(d_sum_weights, n.d_sum_weights);
}

//------------------------------------------------------------------------------------------------//
// end of Norms_Base.cc
//------------------------------------------------------------------------------------------------//
