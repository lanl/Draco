//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   c4/C4_MPI_free_inherited_comm.cc
 * \author Thomas M. Evans
 * \date   Thu Mar 21 16:56:17 2002
 * \brief  C4 MPI implementation.
 * \note   Copyright (C) 2014-2024 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "C4_Functions.hh"
#include "c4/config.h"

#ifdef C4_MPI

namespace rtt_c4 {

//------------------------------------------------------------------------------------------------//
// MPI COMMUNICATOR
//------------------------------------------------------------------------------------------------//

// NOLINTBEGIN(misc-include-cleaner)

extern MPI_Comm communicator;

//------------------------------------------------------------------------------------------------//

void free_inherited_comm() {
  if (communicator != MPI_COMM_WORLD) {
    MPI_Comm_free(&communicator);
    communicator = MPI_COMM_WORLD;
  }
  return;
}

// NOLINTEND(misc-include-cleaner)

} // end namespace rtt_c4

#endif // C4_MPI

//------------------------------------------------------------------------------------------------//
// end of C4_MPI_free_inherited_comm.cc
//------------------------------------------------------------------------------------------------//
