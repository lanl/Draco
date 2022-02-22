//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   Norms_Index.t.cc
 * \author Rob Lowrie
 * \date   Fri Jan 14 13:00:47 2005
 * \brief  Instantiates Norms_Index for some types.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "Norms_Index.t.hh"
#include "Index_Labeled.hh"
#include "Index_Proc.hh"

namespace rtt_norms {

template class Norms_Index<size_t>;
template class Norms_Index<Index_Labeled>;
template class Norms_Index<Index_Proc>;

} // namespace rtt_norms

//------------------------------------------------------------------------------------------------//
// end of norms/Norms_Index.t.cc
//------------------------------------------------------------------------------------------------//
