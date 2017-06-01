//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   c4/C4_Status.hh
 * \author Robert B. Lowrie
 * \date   Friday May 19 6:54:21 2017
 * \brief  C4_Status class definition.
 * \note   Copyright (C) 2016-2017 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//

#ifndef c4_C4_Status_hh
#define c4_C4_Status_hh

// C4 package configure
#include "c4/config.h"
#include "ds++/Assert.hh"

#ifdef C4_MPI
#include "c4_mpi.h"
#endif

namespace rtt_c4 {

//===========================================================================//
/*!
 * \class C4_Status
 * \brief Status container for communications.
 *
 * This class contains the status information for communications.  For MPI,
 * this class wraps MPI_Status.
 */
//===========================================================================//

class DLL_PUBLIC_c4 C4_Status {

#ifdef C4_MPI
  typedef MPI_Status status_type;
#else
  typedef int status_type;
#endif

  status_type d_status;

public:
  // Use default ctor, dtor, assignment

  //! Returns the message size (in bytes) of the last communication.
  int get_message_size() const;

  //! Returns the sending rank of the last communication.
  int get_source() const;

  //! Return a handle to the underlying data status object.
  status_type *get_status_obj() { return &d_status; }
};

} // end namespace rtt_c4

#endif // c4_C4_Status_hh

//---------------------------------------------------------------------------//
// end of c4/C4_Status.hh
//---------------------------------------------------------------------------//
