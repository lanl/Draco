//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   c4/C4_Status.hh
 * \author Robert B. Lowrie
 * \date   Friday May 19 6:54:21 2017
 * \brief  C4_Status class definition.
 * \note   Copyright (C) 2017-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef c4_C4_Status_hh
#define c4_C4_Status_hh

// C4 package configure
#include "c4/config.h"
#include "ds++/Assert.hh"

#ifdef C4_MPI
#include "c4_mpi.h"
#endif

namespace rtt_c4 {

//================================================================================================//
/*!
 * \class C4_Status
 * \brief Status container for communications.
 *
 * This class contains the status information for communications.  For MPI, this class wraps
 * MPI_Status.
 */
//================================================================================================//

class C4_Status {

#ifdef C4_MPI
  using status_type = MPI_Status;
#else
  using status_type = int;
#endif

  status_type d_status;

public:
  //! Returns the message size (in bytes) of the last communication.
  int get_message_size() const;

  //! Returns the sending rank of the last communication.
  int get_source() const;

  //! Return a handle to the underlying data status object.
  status_type *get_status_obj() { return &d_status; }
};

} // end namespace rtt_c4

#endif // c4_C4_Status_hh

//------------------------------------------------------------------------------------------------//
// end of c4/C4_Status.hh
//------------------------------------------------------------------------------------------------//
