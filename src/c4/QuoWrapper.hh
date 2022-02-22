//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   c4/QuoWrapper.hh
 * \author Kelly Thompson
 * \date   Friday, Nov 29, 2019, 18:56 pm
 * \brief  A singleton, light-weight wrapper around Quo (github.com/lanl/libquo)
 * \note   Copyright (C) 2019-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_c4_quowrapper_hh
#define rtt_c4_quowrapper_hh

#include "c4/config.h"
#include "ds++/Assert.hh"
#include <cstdint>
#include <string>

namespace rtt_c4 {

//================================================================================================//
/*!
 * \class QuoWrapper
 * \brief A singlton that is a light-weight wrapper around the Quo API
 *        (https://gitlab.lanl.gov/lanl/libquo).
 */
//================================================================================================//

class QuoWrapper {
public:
  /*! \brief Destructor
   *
   * Because Quo must be torn-down before MPI_finalize is called and because this class uses the
   * RAII pattern, we destroy Quo via QUO_free in C4_MPI.cc's finalize function.  This free's the
   * quo context.
   *
   * This destructor body is only called when the program ends.
   */
  ~QuoWrapper() { quo_free(); }

  /* \brief Destructor helper
   *
   * This function will free the quo context.  It is called by C4_MPI.cc's finalize function.
   *
   * \sa Notes associated with the QuoWrapper destructor.
   */
  static void quo_free() {
#ifdef HAVE_LIBQUO
    if (quo != nullptr)
      QUO_free(quo);
#endif
    return;
  }

  //! disable copy construction
  QuoWrapper(const QuoWrapper &other) = delete;

  //! disable assignment
  QuoWrapper operator=(const QuoWrapper &rhs) = delete;

  //! disable move construction
  QuoWrapper(const QuoWrapper &&other) = delete;

  //! disable move-assignment
  QuoWrapper operator=(const QuoWrapper &&rhs) = delete;

  //! Conruct the singleton on first use.
  static bool is_initialized();

  //! Number of cores per node used by the current process.
  static uint32_t num_cores();

  //! Number of hardware threads per node
  static uint32_t num_hw_threads();

  //! Number of nodes used by the current process.
  static uint32_t num_nodes();

  //! Number of cores per node used by the current process.
  static uint32_t num_mpi_ranks_per_node();

  //! Number of sockets per node
  static uint32_t num_sockets_per_node();

  //! Number of numa domains per node
  static uint32_t num_numanodes_per_node();

  //! Is this process bound to a resource?
  static bool is_bound();

  //! Return a string representation of the caller's current binding policy (cpuset) in a
  //! hexadecimal format.
  static std::string bindings();

  //! Return a pointer (handle) to the static-global QUO_context
#ifdef HAVE_LIBQUO
  static QUO_t *fetch_quo() {
    // Construct on first use.
    Remember(bool isinit =) QuoWrapper::is_initialized();
    Check(isinit);
    return quo;
  }
#endif

private:
  // >> DATA <<

#ifdef HAVE_LIBQUO
  //! Handle to the quo context object
  DLL_PUBLIC_c4 static QUO_context quo;
#endif

  //! Private pointer to this object (this defines the singleton).
  DLL_PUBLIC_c4 static QuoWrapper *instance;

  // >> MANIPULATORS <<

  //! Private constructor so that no objects can be created.
  QuoWrapper();

}; // end class QuoWrapper

} // namespace rtt_c4

#endif // rtt_c4_quowrapper_hh

//------------------------------------------------------------------------------------------------//
// end of QuoWrapper.hh
//------------------------------------------------------------------------------------------------//
