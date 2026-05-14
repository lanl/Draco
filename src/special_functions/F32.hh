//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   special_functions/F32.hh
 * \author Kent Budge
 * \brief  Compute Fermi-Dirac function of 1/2 order
 * \note   Copyright (C) 2010-2026 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef sf_F32_hh
#define sf_F32_hh

#include "dsxx/config.h"

namespace rtt_sf {
//! Calculate Fermi-Dirac integral of index 3/2.
template <typename OrderedField> OrderedField F32(OrderedField const &x);

} // end namespace rtt_sf

#endif // sf_F32_hh

//------------------------------------------------------------------------------------------------//
// end of sf/F32.hh
//------------------------------------------------------------------------------------------------//
