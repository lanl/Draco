//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   Norms_Index_pt.cc
 * \author Rob Lowrie
 * \date   Fri Jan 14 13:00:47 2005
 * \brief  Instantiates Norms_Index for some types.
 * \note   Copyright (C) 2024 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "norms/Index_Labeled.hh"
#include "norms/Index_Proc.hh"
#include "norms/Norms_Index.hh"
#include "norms/Norms_Index.t.hh" // NOLINT(misc-include-cleaner)

namespace rtt_norms {

template class Norms_Index<size_t>;
template class Norms_Index<Index_Labeled>;
template class Norms_Index<Index_Proc>;

} // namespace rtt_norms

//------------------------------------------------------------------------------------------------//
// end of norms/Norms_Index.t.cc
//------------------------------------------------------------------------------------------------//
