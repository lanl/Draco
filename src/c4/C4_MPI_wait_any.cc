//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   c4/C4_MPI_wait_any.cc
 * \author Thomas M. Evans
 * \date   Thu Mar 21 16:56:17 2002
 * \brief  C4 MPI implementation.
 * \note   Copyright (C) 2014-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "c4/config.h"
#include <vector>

#ifdef C4_MPI

#include "C4_Req.hh"

namespace rtt_c4 {

//------------------------------------------------------------------------------------------------//
unsigned wait_any(unsigned count, C4_Req *requests) {
  using std::vector;

  vector<MPI_Request> array_of_requests(count);
  for (unsigned i = 0; i < count; ++i) {
    if (requests[i].inuse())
      array_of_requests[i] = requests[i].r();
    else
      array_of_requests[i] = MPI_REQUEST_NULL;
  }
  int index;
  MPI_Waitany(count, &array_of_requests[0], &index, MPI_STATUSES_IGNORE);
  requests[index] = C4_Req();

  return index;
}

} // end namespace rtt_c4

#endif // C4_MPI

//------------------------------------------------------------------------------------------------//
// end of C4_MPI_wait_any.cc
//------------------------------------------------------------------------------------------------//
