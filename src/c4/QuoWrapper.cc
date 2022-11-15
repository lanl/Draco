//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   c4/QuoWrapper.cc
 * \author Kelly Thompson
 * \date   Friday, Nov 29, 2019, 18:15 pm
 * \brief  Wrap libquo (github.com/lanl/libquo)
 * \note   Copyright (C) 2019-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "QuoWrapper.hh"
#include "C4_Functions.hh"
#include "ds++/Assert.hh"

namespace rtt_c4 {

//------------------------------------------------------------------------------------------------//
/*! \brief Initialize pointer to zero so that it can be initialized in first call to getInstance
 */
QuoWrapper *QuoWrapper::instance = nullptr;
#ifdef HAVE_LIBQUO
QUO_context QuoWrapper::quo = nullptr;
#endif

//------------------------------------------------------------------------------------------------//
/*!
 * \brief default constructor (private!)
 *
 * Use a RAII model - construction and initialization of the quo objects and the data for this
 * wrapper.
 *
 * Singleton: This function is private and can only be called by the static is_initialized member
 * function.
 */
QuoWrapper::QuoWrapper() {
  Insist(isMpiInit(), "Cannot create a QuoWrapper unless MPI is already initialized");
#ifdef HAVE_LIBQUO
  if (quo == nullptr) {
    Remember(auto retval =) QUO_create(&quo, MPI_COMM_WORLD);
    Check(retval == QUO_SUCCESS);
  }
#endif
}

//------------------------------------------------------------------------------------------------//
/*! \brief Calling this public function should always return true.
 *
 * If this singleton has not be created, then it will be created. Otherwise, just check that the
 * pointer to instance is valid and return true.
 *
 * \pre Must be called (constructed) after MPI_Init!
 */
bool QuoWrapper::is_initialized() {
  if (instance == nullptr)
    instance = new QuoWrapper;
  return instance != nullptr;
}

//------------------------------------------------------------------------------------------------//
//! Number of cores per node
uint32_t QuoWrapper::num_cores() {
  // Construct on first use.
  Remember(bool isinit =) QuoWrapper::is_initialized();
  Check(isinit);
  int n(-1);
#ifdef HAVE_LIBQUO
  Remember(auto retval =) QUO_ncores(quo, &n);
  Check(QUO_SUCCESS == retval);
  Ensure(n > 0);
#else
  n = 0;
#endif
  return static_cast<uint32_t>(n);
}

//! Number of hardware threads
uint32_t QuoWrapper::num_hw_threads() {
  // Construct on first use.
  Remember(bool isinit =) QuoWrapper::is_initialized();
  Check(isinit);
  int n(-1);
#ifdef HAVE_LIBQUO
  Remember(auto retval =) QUO_npus(quo, &n);
  Check(QUO_SUCCESS == retval);
  Ensure(n > 0);
#else
  n = 0;
#endif
  return static_cast<uint32_t>(n);
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Number of nodes used by the current process
 *
 * - A 'node' is a machine that is connected to others via MPI or high speed network. Each node has
 *   a number of sockets, numa nodes, and a fixed number of cores.
 * - Report how many of these nodes are used by the current job.  If our allocation has 2 nodes, but
 *   all MPI ranks are assigned to the same node, this will report 1.
 */
uint32_t QuoWrapper::num_nodes() {
  // Construct on first use.
  Remember(bool isinit =) QuoWrapper::is_initialized();
  Check(isinit);
  int n(-1);
#ifdef HAVE_LIBQUO
  Remember(auto retval =) QUO_nnodes(quo, &n);
  Check(QUO_SUCCESS == retval);
  Ensure(n > 0);
#else
  n = 0;
#endif
  return static_cast<uint32_t>(n);
}

//------------------------------------------------------------------------------------------------//
//! Number of cores per node used by the current process.
uint32_t QuoWrapper::num_mpi_ranks_per_node() {
  // Construct on first use.
  Remember(bool isinit =) QuoWrapper::is_initialized();
  Check(isinit);
  int n(-1);
#ifdef HAVE_LIBQUO
  Remember(auto retval =) QUO_nqids(quo, &n);
  Check(QUO_SUCCESS == retval);
  Ensure(n > 0);
#else
  n = 0;
#endif
  return static_cast<uint32_t>(n);
}

//------------------------------------------------------------------------------------------------//
//! Number of sockets per node
uint32_t QuoWrapper::num_sockets_per_node() {
  // Construct on first use.
  Remember(bool isinit =) QuoWrapper::is_initialized();
  Check(isinit);
  int n(-1);
#ifdef HAVE_LIBQUO
  Remember(auto retval =) QUO_nsockets(quo, &n);
  Check(QUO_SUCCESS == retval);
#endif
  Ensure(n > 0);
  return static_cast<uint32_t>(n);
}

//------------------------------------------------------------------------------------------------//
//! Number of numa domains per node
uint32_t QuoWrapper::num_numanodes_per_node() {
  // Construct on first use.
  Remember(bool isinit =) QuoWrapper::is_initialized();
  Check(isinit);
  int n(-1);
#ifdef HAVE_LIBQUO
  Remember(auto retval =) QUO_nnumanodes(quo, &n);
  Check(QUO_SUCCESS == retval);
  // If n is set to zero above and no error is returned, then set n to 1.  This seems to occur for
  // WLS2.
  if (n == 0)
    n = 1;
  Ensure(n > 0);
#else
  n = 0;
#endif
  return static_cast<uint32_t>(n);
}

//------------------------------------------------------------------------------------------------//
//! Is this process bound to a resource?
bool QuoWrapper::is_bound() {
  // Construct on first use.
  Remember(bool isinit =) QuoWrapper::is_initialized();
  Check(isinit);
  int n(-1);
#ifdef HAVE_LIBQUO
  Remember(auto retval =) QUO_bound(quo, &n);
  Check(QUO_SUCCESS == retval);
  Check(n == 0 || n == 1);
#endif
  return (n == 1);
}

//------------------------------------------------------------------------------------------------//
//! Return a string representation of the caller's current binding policy ! (cpuset) in a
//hexadecimal format.
std::string QuoWrapper::bindings() {
  // Construct on first use.
  Remember(bool isinit =) QuoWrapper::is_initialized();
  Check(isinit);
#ifdef HAVE_LIBQUO
  char *cbindstr = nullptr;
  Remember(auto retval =) QUO_stringify_cbind(quo, &cbindstr);
  Check(QUO_SUCCESS == retval);
  std::string cppbindstr(cbindstr);
  free(cbindstr); // NOLINT [hicpp-no-malloc]
#else
  std::string cppbindstr("unknown - quo not available");
#endif
  return cppbindstr;
}

} // end namespace rtt_c4

//------------------------------------------------------------------------------------------------//
// end of QuoWrapper.cc
//------------------------------------------------------------------------------------------------//
