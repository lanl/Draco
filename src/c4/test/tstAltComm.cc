//--------------------------------------------*-C++-*---------------------------------------------//
/*!
  * \file   c4/test/tstAltComm.cc
  * \author Mathew Cleveland 
  * \brief  Test alternative (non MPI_COMM_WORLD) communictor
  * \note   Copyright (C) 2025-2026 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "dsxx/Assert.hh"
#include "dsxx/Release.hh"
#include "dsxx/UnitTest.hh"
#include "c4/C4_Functions.hh"
#include "c4/C4_MPI.hh"
#include "c4/ParallelUnitTest.hh"
#include <iostream>
#include <mpi.h>
#include <sstream>

namespace {

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//
using namespace std;

MPI_Comm set_alt_comm(rtt_dsxx::UnitTest &ut, bool &comm_rank) {

  FAIL_IF_NOT(rtt_c4::nodes() == 4);

  // Store for checking
  const int original_node_id{rtt_c4::node()};

  // Set the communication color=0 non-comm color=1 comm.
  int color{0};
  if (rtt_c4::node() > 1)
    color = 1;
  MPI_Comm old_comm{rtt_c4::communicator};
  MPI_Comm new_comm;
  MPI_Comm_split(rtt_c4::communicator, color, original_node_id, &new_comm);
  rtt_c4::communicator = new_comm;

  // These calls would hang if they were called on color==0 ranks
  comm_rank = false;
  if (color == 1) {
    comm_rank = true;
    rtt_c4::global_barrier();
    FAIL_IF_NOT(rtt_c4::nodes() == 2);
    FAIL_IF_NOT((original_node_id == 2 && rtt_c4::node() == 0) ||
                (original_node_id == 3 && rtt_c4::node() == 1));
    int val{1};
    rtt_c4::global_sum(val);
    FAIL_IF_NOT(val == 2);
  }

  if (ut.numFails == 0) {
    ostringstream m;
    m << "Set alternative communicator passed on " << rtt_c4::node();
    PASSMSG(m.str());
  } else {
    ostringstream m;
    m << "set alternative communicator failed on " << rtt_c4::node();
    FAILMSG(m.str());
  }
  return old_comm;
}

void test_persistent_comm(rtt_dsxx::UnitTest &ut, bool &comm_rank) {
  if (!comm_rank)
    return;
  rtt_c4::global_barrier();
  FAIL_IF_NOT(rtt_c4::nodes() == 2);
  int val{1};
  rtt_c4::global_sum(val);
  FAIL_IF_NOT(val == 2);

  if (ut.numFails == 0) {
    ostringstream m;
    m << " Persistent alternative communicator passed on " << rtt_c4::node();
    PASSMSG(m.str());
  } else {
    ostringstream m;
    m << " Persistent alternative communicator failed on " << rtt_c4::node();
    FAILMSG(m.str());
  }
  rtt_c4::global_barrier();
  return;
}
void test_world_comm(rtt_dsxx::UnitTest &ut) {
  rtt_c4::global_barrier();
  FAIL_IF_NOT(rtt_c4::nodes() == 4);
  int val{1};
  rtt_c4::global_sum(val);
  FAIL_IF_NOT(val == 4);

  if (ut.numFails == 0) {
    ostringstream m;
    m << " World communicator passed on " << rtt_c4::node();
    PASSMSG(m.str());
  } else {
    ostringstream m;
    m << " World communicator failed on " << rtt_c4::node();
    FAILMSG(m.str());
  }
  rtt_c4::global_barrier();
  return;
}

} // namespace

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  rtt_c4::ParallelUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    std::cout << "This is " << rtt_c4::get_processor_name() << '\n';
    Insist(rtt_c4::nodes() == 4, "This test is designed for exactly 4 processors.");
    test_world_comm(ut);
    bool comm_rank(true);
    MPI_Comm old_comm = set_alt_comm(ut, comm_rank);
    test_persistent_comm(ut, comm_rank);
    // Free Alternative communicator and reset original communicator and test world comm again.
    MPI_Comm_free(&rtt_c4::communicator);
    rtt_c4::communicator = old_comm;
    test_world_comm(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstAltComm.cc
//------------------------------------------------------------------------------------------------//
