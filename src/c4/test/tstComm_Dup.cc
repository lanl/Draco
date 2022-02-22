//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   c4/test/tstComm_Dup.cc
 * \author Thomas M. Evans
 * \date   Thu Jul 18 11:10:10 2002
 * \brief  test Communicator Duplication
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

#include "c4/ParallelUnitTest.hh"
#include "ds++/Release.hh"
#include <sstream>

using namespace std;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

void test_mpi_comm_dup(rtt_dsxx::UnitTest &ut) {

// we only run this particular test when MPI is on
#ifdef C4_MPI

  Require(rtt_c4::nodes() == 4);

  int node = rtt_c4::node();
  int snode = 0;

  // split up nodes (two communicators) 0 -> 0, 2 -> 1 and 1 -> 0, 3 -> 1
#if defined(OMPI_MAJOR_VERSION)
  MPI_Comm new_comm = nullptr;
#else
  // if defined(MSMPI_VER) || defined(CRAY_MPICH_VERSION)
  MPI_Comm new_comm = 0; // MS_MPI: 'typedef int MPI_Comm;'
#endif

  if (node == 1) {
    MPI_Comm_split(MPI_COMM_WORLD, 0, 0, &new_comm);
    MPI_Comm_rank(new_comm, &snode);
  } else if (node == 3) {
    MPI_Comm_split(MPI_COMM_WORLD, 0, 0, &new_comm);
    MPI_Comm_rank(new_comm, &snode);
  } else {
    MPI_Comm_split(MPI_COMM_WORLD, 1, 0, &new_comm);
    MPI_Comm_rank(new_comm, &snode);
  }

  // we haven't set the communicator yet so we should still have 4 nodes
  FAIL_IF_NOT(rtt_c4::nodes() == 4);

  // now duplicate the communicator on each processor
  rtt_c4::inherit(new_comm);

  // each processor should see two nodes
  FAIL_IF_NOT(rtt_c4::nodes() == 2);

  // test data send/receive
  int data = 0;

  // do some tests on each processor
  if (node == 0) {

    FAIL_IF_NOT(rtt_c4::node() == 0);
    // set data to 10 and send it out
    data = 10;
    rtt_c4::send(&data, 1, 1, 100);

  } else if (node == 1) {

    FAIL_IF_NOT(rtt_c4::node() == 0);
    // set data to 20 and send it out
    data = 20;
    rtt_c4::send(&data, 1, 1, 100);

  } else if (node == 2) {

    FAIL_IF_NOT(rtt_c4::node() == 1);
    FAIL_IF_NOT(data == 0);
    rtt_c4::receive(&data, 1, 0, 100);
    FAIL_IF_NOT(data == 10);

  } else if (node == 3) {

    FAIL_IF_NOT(rtt_c4::node() == 1);
    FAIL_IF_NOT(data == 0);
    rtt_c4::receive(&data, 1, 0, 100);
    FAIL_IF_NOT(data == 20);
  }

  // now free the communicator on each processor
  rtt_c4::free_inherited_comm();

  // the free should have set back to COMM_WORLD
  FAIL_IF_NOT(rtt_c4::nodes() == 4);

  rtt_c4::global_barrier();
  if (ut.numFails == 0) {
    std::ostringstream m;
    m << "Communicator duplicated successfully on " << rtt_c4::node();
    PASSMSG(m.str());
  }

  rtt_c4::global_barrier();
  MPI_Comm_free(&new_comm);

#else
  /* DRACO_C4 = SCALAR */
  PASSMSG("Nothing to do since DRACO_C4 = SCALAR");

#endif

  return;
}

//------------------------------------------------------------------------------------------------//
void test_comm_dup(rtt_dsxx::UnitTest &ut) {
// we only run this test scalar
#ifdef C4_SCALAR

  int node = rtt_c4::node();

  // now duplicate the communicator on each processor
  rtt_c4::inherit(node);

  // check the number of nodes
  FAIL_IF_NOT(rtt_c4::nodes() == 1);

  rtt_c4::free_inherited_comm();

  // check the number of nodes
  FAIL_IF_NOT(rtt_c4::nodes() == 1);

  if (ut.numFails == 0)
    PASSMSG("Scalar Comm duplication/free works okay.");

#endif

// check duping/freeing MPI_COMM_WORLD
#ifdef C4_MPI

  int nodes = rtt_c4::nodes();

  MPI_Comm comm_world = MPI_COMM_WORLD;
  rtt_c4::inherit(comm_world);

  FAIL_IF_NOT(rtt_c4::nodes() == nodes);

  // try a global sum to check
  int x = 10;
  rtt_c4::global_sum(x);
  FAIL_IF_NOT(x == 10 * nodes);

  rtt_c4::free_inherited_comm();

  // we should be back to comm world
  FAIL_IF_NOT(rtt_c4::nodes() == nodes);

  // try a global sum to check
  int y = 20;
  rtt_c4::global_sum(y);
  FAIL_IF_NOT(y == 20 * nodes);

  if (ut.numFails == 0)
    PASSMSG("MPI_COMM_WORLD Comm duplication/free works okay.");

#endif
  return;
}

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  rtt_c4::ParallelUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    // >>> UNIT TESTS
    if (rtt_c4::nodes() == 4)
      test_mpi_comm_dup(ut);

    test_comm_dup(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstComm_Dup.cc
//------------------------------------------------------------------------------------------------//
