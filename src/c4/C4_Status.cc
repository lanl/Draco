//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   c4/C4_Status.cc
 * \author Robert B. Lowrie
 * \date   Friday May 19 9:31:33 2017
 * \brief  C4_Status member definitions.
 * \note   Copyright (C) 2017-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "C4_Status.hh"

namespace rtt_c4 {

//------------------------------------------------------------------------------------------------//
int C4_Status::get_message_size() const {
  int ms(0); // return value
#ifdef C4_MPI
  MPI_Get_count(&d_status, MPI_CHAR, &ms);
#endif
  return ms;
}

//------------------------------------------------------------------------------------------------//
int C4_Status::get_source() const {
  int s(0); // return value
#ifdef C4_MPI
  s = d_status.MPI_SOURCE;
#endif
  return s;
}

} // end namespace rtt_c4

//------------------------------------------------------------------------------------------------//
// end of C4_Status.cc
//------------------------------------------------------------------------------------------------//
