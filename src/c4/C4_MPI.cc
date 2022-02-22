//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   c4/C4_MPI.cc
 * \author Thomas M. Evans
 * \date   Thu Mar 21 16:56:17 2002
 * \brief  C4 MPI implementation.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "c4/config.h"
#include <vector>

#ifdef C4_MPI

#include "C4_Functions.hh"
#include "C4_Req.hh"
#include "C4_sys_times.h"
#include "QuoWrapper.hh"

namespace rtt_c4 {

//------------------------------------------------------------------------------------------------//
// MPI COMMUNICATOR
//------------------------------------------------------------------------------------------------//

MPI_Comm communicator = MPI_COMM_WORLD;
bool initialized(false);

//------------------------------------------------------------------------------------------------//
// Any source rank
//------------------------------------------------------------------------------------------------//

const int any_source = MPI_ANY_SOURCE;

//------------------------------------------------------------------------------------------------//
// Null source/destination rank
//------------------------------------------------------------------------------------------------//

const int proc_null = MPI_PROC_NULL;

//------------------------------------------------------------------------------------------------//
// SETUP FUNCTIONS
//------------------------------------------------------------------------------------------------//

int initialize(int &argc, char **&argv, int required) {
  int provided;
  int result = MPI_Init_thread(&argc, &argv, required, &provided);
  initialized = (result == MPI_SUCCESS);
  Check(initialized);

  // Resync clocks for Darwin mpich
  Remember(double foo(MPI_Wtick()););
  Ensure(foo > 0.0);
  return provided;
}

//------------------------------------------------------------------------------------------------//

void finalize() {
  // If Libquo is active, it must be torn-down before MPI_Finalize is called.  Otherwise, this call
  // is a no-op.
  QuoWrapper::quo_free();
  MPI_Finalize();
  return;
}

//------------------------------------------------------------------------------------------------//

void type_free(C4_Datatype &old_type) { MPI_Type_free(&old_type); }

//------------------------------------------------------------------------------------------------//
// QUERY FUNCTIONS
//------------------------------------------------------------------------------------------------//

int node() {
  int node = 1;
  MPI_Comm_rank(communicator, &node);
  Check(node >= 0);
  return node;
}
uint32_t rank() {
  int rank = 1;
  MPI_Comm_rank(communicator, &rank);
  Check(rank >= 0);
  return static_cast<uint32_t>(rank);
}

//------------------------------------------------------------------------------------------------//

int nodes() {
  int nodes = 0;
  MPI_Comm_size(communicator, &nodes);
  Check(nodes > 0);
  return nodes;
}
uint32_t nranks() {
  int nranks = 0;
  MPI_Comm_size(communicator, &nranks);
  Check(nranks > 0);
  return static_cast<uint32_t>(nranks);
}

//------------------------------------------------------------------------------------------------//
// BARRIER FUNCTIONS
//------------------------------------------------------------------------------------------------//

void global_barrier() {
  MPI_Barrier(communicator);
  return;
}

//------------------------------------------------------------------------------------------------//
// TIMING FUNCTIONS
//------------------------------------------------------------------------------------------------//

// overloaded function (no arguments)
double wall_clock_time() { return MPI_Wtime(); }
// overloaded function (provide POSIX timer information).
double wall_clock_time(DRACO_TIME_TYPE &now) {
// obtain POSIX timer information and return it to the user via the reference value argument "now".
#ifdef _MSC_VER
  now = std::chrono::high_resolution_clock::now();
#else
  times(&now);
#endif
  // This function will return the MPI wall-clock time.
  return MPI_Wtime();
}

//------------------------------------------------------------------------------------------------//
double wall_clock_resolution() { return MPI_Wtick(); }

//------------------------------------------------------------------------------------------------//
// PROBE/WAIT FUNCTIONS
//------------------------------------------------------------------------------------------------//

bool probe(int source, int tag, int &message_size) {
  //! \todo Change message_size to C4_Status to allow source = any_source
  //Require(source == any_source || (source >= 0 && source < nodes()));
  Require(source >= 0 && source < nodes());

  int flag;
  MPI_Status status;

  // post a non-blocking probe
  MPI_Iprobe(source, tag, communicator, &flag, &status);

  if (!flag)
    return false;

  MPI_Get_count(&status, MPI_CHAR, &message_size);

  return true;
}

//------------------------------------------------------------------------------------------------//
void blocking_probe(int source, int tag, int &message_size) {
  //! \todo Change message_size to C4_Status to allow source = any_source
  //Require(source == any_source || (source >= 0 && source < nodes()));
  Require(source >= 0 && source < nodes());

  MPI_Status status;
  MPI_Probe(source, tag, communicator, &status);
  MPI_Get_count(&status, MPI_CHAR, &message_size);
}

//------------------------------------------------------------------------------------------------//
void wait_all(const unsigned count, C4_Req *const requests) {

  // Nothing to do if count is zero.
  if (count == 0)
    return;

  std::vector<MPI_Request> array_of_requests(count);
  for (unsigned i = 0; i < count; ++i) {
    if (requests[i].inuse())
      array_of_requests[i] = requests[i].r();
    else
      array_of_requests[i] = MPI_REQUEST_NULL;
  }
  Remember(int check =) MPI_Waitall(count, &array_of_requests[0], MPI_STATUSES_IGNORE);
  Check(check == MPI_SUCCESS);
  return;
}

//------------------------------------------------------------------------------------------------//
std::vector<int> wait_all_with_source(const unsigned count, C4_Req *const requests) {

  // Nothing to do if count is zero.
  if (count == 0)
    return {};

  // Return value -- rank IDs for all message sources.
  std::vector<int> msg_sources(count);

  std::vector<MPI_Request> array_of_requests(count);
  // This is needed to obtain the source rank for each message:
  std::vector<MPI_Status> array_of_statuses(count);
  for (unsigned i = 0; i < count; ++i) {
    if (requests[i].inuse())
      array_of_requests[i] = requests[i].r();
    else
      array_of_requests[i] = MPI_REQUEST_NULL;
  }
  Remember(int check =) MPI_Waitall(count, &array_of_requests[0], &array_of_statuses[0]);
  Check(check == MPI_SUCCESS);

  for (unsigned p = 0; p < count; p++) {
    msg_sources[p] = array_of_statuses[p].MPI_SOURCE;
  }

  return msg_sources;
}

//------------------------------------------------------------------------------------------------//
// ABORT
//------------------------------------------------------------------------------------------------//
int abort(int error) {
  int rerror = MPI_Abort(communicator, error);
  return rerror;
}

//------------------------------------------------------------------------------------------------//
// Helpers
//------------------------------------------------------------------------------------------------//
bool isScalar() { return !initialized; }
bool isMpiInit() { return initialized; }

} // end namespace rtt_c4

//! Set c4's initialzed variable to true (called from Fortran tests)
void setMpiInit() {
  rtt_c4::initialized = true;
  return;
}

#endif // C4_MPI

//------------------------------------------------------------------------------------------------//
// end of C4_MPI.cc
//------------------------------------------------------------------------------------------------//
