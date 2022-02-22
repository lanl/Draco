//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   rng/Sobol_Sequence.cc
 * \author Kent Budge
 * \date   Thu Dec 22 13:38:35 2005
 * \brief  Implementation file for Sobol Sequence.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "Sobol_Sequence.hh"

namespace rtt_rng {

//------------------------------------------------------------------------------------------------//
/*!
 * \param dimension Dimension of the subrandom vector returned by this object.
 */
Sobol_Sequence::Sobol_Sequence(unsigned const dimension)
    : gsl_(gsl_qrng_alloc(gsl_qrng_sobol, dimension)), values_(dimension) {
  Require(dimension > 0);
  shift();
  Ensure(check_class_invariants());
}

//------------------------------------------------------------------------------------------------//
void Sobol_Sequence::shift() {
  gsl_qrng_get(gsl_, &values_[0]);
  Ensure(check_class_invariants());
}

//------------------------------------------------------------------------------------------------//
Sobol_Sequence::~Sobol_Sequence() { gsl_qrng_free(gsl_); }

//------------------------------------------------------------------------------------------------//
bool Sobol_Sequence::check_class_invariants() const { return gsl_ != nullptr; }

} // end namespace rtt_rng

//------------------------------------------------------------------------------------------------//
// end of Sobol_Sequence.cc
//------------------------------------------------------------------------------------------------//
