//--------------------------------------------*-C++-*---------------------------------------------//
/*!
  \file   Norms.hh
  \author Rob Lowrie
  \date   Fri Jan 14 13:00:32 2005
  \brief  Header file for Norms.
  \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_norms_Norms_hh
#define rtt_norms_Norms_hh

#include "Norms_Index.hh"

namespace rtt_norms {

// Norms is the most common usage of Norms_Index.  A default template argument
// of size_t for Norms_Indexis not used, because of the horrible syntax
// Norms_Index<>.

//! Convenience definition for Norms_Index<size_t>.
using Norms = Norms_Index<size_t>;

} // namespace rtt_norms

#endif // rtt_norms_Norms_hh

//------------------------------------------------------------------------------------------------//
// end of norms/Norms.hh
//------------------------------------------------------------------------------------------------//
