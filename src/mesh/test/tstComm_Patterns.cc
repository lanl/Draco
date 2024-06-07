//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   mesh/test/tstComm_Patterns.cc
 * \author Mathew Cleveland
 * \date   Nov. 2023
 * \brief  Comm_Patterns test.
 * \note   Copyright (C) 2023-2024 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "Test_Topology.hh"
#include "dsxx/Release.hh"
#include "c4/ParallelUnitTest.hh"
#include "mesh/Comm_Patterns.hh"
#include "mesh/DD_Topology.hh"
#include "mesh/Rep_Topology.hh"
#include <sstream>

using namespace std;

using rtt_mesh::Comm_Patterns;
using rtt_mesh::DD_Topology;
using rtt_mesh::Rep_Topology;
using rtt_mesh::Topology;

//------------------------------------------------------------------------------------------------//
// TEST HELPERS
// ------------------------------------------------------------------------------------------------//

// Build a 2 processor 8 cell topology illustrated in autodoc/DD_Comm_Patterns.fig.

std::shared_ptr<Topology> build_2p_DD_Topology() {
  Require(rtt_c4::nodes() == 2);

  std::shared_ptr<Topology> topology;

  // num_global = 8;
  vector<int> bc;
  vector<int> bp;

  vector<int> cpp(2, 0);
  cpp[0] = 3;
  cpp[1] = 5;

  int n = rtt_c4::node();

  if (n == 0) {
    bc.resize(4);
    bp.resize(4, 1);

    bc[0] = 2;
    bc[1] = 3;
    bc[2] = 5;
    bc[3] = 4;
  } else if (n == 1) {
    bc.resize(3);
    bp.resize(3, 0);

    bc[0] = 2;
    bc[1] = 3;
    bc[2] = 1;
  }

  topology = std::make_shared<DD_Topology>(cpp, bc, bp);
  return topology;
}

//------------------------------------------------------------------------------------------------//
// Build a 4 process 8 cell decomposition with multiple boundary cell indices.  See
// autodoc/DD_Comm_Patterns_multiple.fig.

std::shared_ptr<Topology> build_4p_DD_Topology_multiple() {
  Require(rtt_c4::nodes() == 4);

  std::shared_ptr<Topology> topology;

  // num_global = 8;
  vector<int> bc;
  vector<int> bp;

  vector<int> cpp(4, 0);
  cpp[0] = 2;
  cpp[1] = 1;
  cpp[2] = 2;
  cpp[3] = 3;

  int n = rtt_c4::node();

  if (n == 0) {
    bc.resize(4);
    bp.resize(4);

    bc[0] = 2;
    bc[1] = 2;
    bc[2] = 1;
    bc[3] = 2;

    bp[0] = 3;
    bp[1] = 2;
    bp[2] = 3;
    bp[3] = 2;
  } else if (n == 1) {
    bc.resize(6);
    bp.resize(6);

    bc[0] = 1;
    bc[1] = 1;
    bc[2] = 2;
    bc[3] = 3;
    bc[4] = 2;
    bc[5] = 1;

    bp[0] = 3;
    bp[1] = 2;
    bp[2] = 2;
    bp[3] = 3;
    bp[4] = 2;
    bp[5] = 2;
  } else if (n == 2) {
    bc.resize(8);
    bp.resize(8);

    bc[0] = 3;
    bc[1] = 1;
    bc[2] = 2;
    bc[3] = 1;
    bc[4] = 2;
    bc[5] = 1;
    bc[6] = 1;
    bc[7] = 1;

    bp[0] = 3;
    bp[1] = 1;
    bp[2] = 0;
    bp[3] = 3;
    bp[4] = 0;
    bp[5] = 3;
    bp[6] = 3;
    bp[7] = 0;
  } else if (n == 3) {
    bc.resize(4);
    bp.resize(4);

    bc[0] = 1;
    bc[1] = 2;
    bc[2] = 1;
    bc[3] = 1;

    bp[0] = 2;
    bp[1] = 0;
    bp[2] = 0;
    bp[3] = 1;
  }

  topology = std::make_shared<DD_Topology>(cpp, bc, bp);
  return topology;
}

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//
// Comm_Patterns DD Test on 8 cell mesh in autodoc/DD_Comm_Patterns.fig.

void DD_Comm_Patterns_2(rtt_dsxx::UnitTest &ut) {
  if (rtt_c4::nodes() != 2)
    return;

  int n = rtt_c4::node();
  int p = 0;
  Comm_Patterns::const_iterator itr;
  vector<int> cells;

  // build a DD topology for a nine cell mesh
  std::shared_ptr<Topology> topology = build_2p_DD_Topology();
  FAIL_IF_NOT(topology->get_parallel_scheme() == "DD");

  // make a Comm_Patterns object
  Comm_Patterns cp;

  // nothing should be set yet
  FAIL_IF(cp);

  // calculate comm_patterns
  cp.calc_patterns(topology);
  FAIL_IF_NOT(cp);

  if (n == 0) {
    if (cp.get_num_recv_procs() == 1)
      FAIL_IF_NOT(cp.get_num_send_procs() == 1);

    // receive maps
    itr = cp.get_recv_begin();

    p = itr->first;
    cells = itr->second;

    FAIL_IF_NOT(p == 1);
    FAIL_IF_NOT(cells.size() == 4);
    FAIL_IF_NOT(cells[0] == 1);
    FAIL_IF_NOT(cells[1] == 2);
    FAIL_IF_NOT(cells[2] == 3);
    FAIL_IF_NOT(cells[3] == 4);

    itr++;
    FAIL_IF_NOT(itr == cp.get_recv_end());

    // send maps
    itr = cp.get_send_begin();

    p = itr->first;
    cells = itr->second;

    FAIL_IF_NOT(p == 1);
    FAIL_IF_NOT(cells.size() == 3);
    FAIL_IF_NOT(cells[0] == 2);
    FAIL_IF_NOT(cells[1] == 3);
    FAIL_IF_NOT(cells[2] == 1);

    itr++;
    FAIL_IF_NOT(itr == cp.get_send_end());
  } else if (n == 1) {
    FAIL_IF_NOT(cp.get_num_recv_procs() == 1);
    FAIL_IF_NOT(cp.get_num_send_procs() == 1);

    // receive maps
    itr = cp.get_recv_begin();

    p = itr->first;
    cells = itr->second;

    FAIL_IF_NOT(p == 0);
    FAIL_IF_NOT(cells.size() == 3);
    FAIL_IF_NOT(cells[0] == 1);
    FAIL_IF_NOT(cells[1] == 2);
    FAIL_IF_NOT(cells[2] == 3);

    itr++;
    FAIL_IF_NOT(itr == cp.get_recv_end());

    // send maps
    itr = cp.get_send_begin();

    p = itr->first;
    cells = itr->second;

    FAIL_IF_NOT(p == 0);
    FAIL_IF_NOT(cells.size() == 4);
    FAIL_IF_NOT(cells[0] == 2);
    FAIL_IF_NOT(cells[1] == 3);
    FAIL_IF_NOT(cells[2] == 5);
    FAIL_IF_NOT(cells[3] == 4);

    itr++;
    FAIL_IF_NOT(itr == cp.get_send_end());
  }

  if (ut.numFails == 0) {
    ostringstream m;
    m << "DD Comm_Patterns 2 proc test ok on " << n;
    PASSMSG(m.str());
  }
}
//------------------------------------------------------------------------------------------------//
// Test Comm_Patterns when a boundary cell has multiple indices.  See
// autodoc/DD_Comm_Patterns_multiple.fig.

void DD_Comm_Patterns_4_multiple_indices(rtt_dsxx::UnitTest &ut) {
  if (rtt_c4::nodes() != 4)
    return;

  int n = rtt_c4::node();
  int p = 0;
  Comm_Patterns::const_iterator itr;
  vector<int> cells;

  // build a DD topology for a nine cell mesh
  std::shared_ptr<Topology> topology = build_4p_DD_Topology_multiple();
  if (topology->get_parallel_scheme() != "DD")
    ITFAILS;

  // make a cp object
  Comm_Patterns cp;

  // nothing should be set yet
  if (cp)
    ITFAILS;

  // calculate comm_patterns
  cp.calc_patterns(topology);
  if (!cp)
    ITFAILS;

  if (n == 0) {
    FAIL_IF_NOT(cp.get_num_recv_procs() == 2);
    FAIL_IF_NOT(cp.get_num_send_procs() == 2);

    // receive maps
    itr = cp.get_recv_begin();

    p = itr->first;
    cells = itr->second;

    FAIL_IF_NOT(p == 2);
    FAIL_IF_NOT(cells.size() == 2);
    FAIL_IF_NOT(cells[0] == 2);
    FAIL_IF_NOT(cells[1] == 4);

    itr++;

    p = itr->first;
    cells = itr->second;

    FAIL_IF_NOT(p == 3);
    FAIL_IF_NOT(cells.size() == 2);
    FAIL_IF_NOT(cells[0] == 1);
    FAIL_IF_NOT(cells[1] == 3);

    itr++;
    FAIL_IF_NOT(itr == cp.get_recv_end());

    // send maps
    itr = cp.get_send_begin();

    p = itr->first;
    cells = itr->second;

    FAIL_IF_NOT(p == 2);
    FAIL_IF_NOT(cells.size() == 3);
    FAIL_IF_NOT(cells[0] == 2);
    FAIL_IF_NOT(cells[1] == 2);
    FAIL_IF_NOT(cells[2] == 1);

    itr++;

    p = itr->first;
    cells = itr->second;

    FAIL_IF_NOT(p == 3);
    FAIL_IF_NOT(cells.size() == 2);
    FAIL_IF_NOT(cells[0] == 2);
    FAIL_IF_NOT(cells[1] == 1);

    itr++;
    FAIL_IF_NOT(itr == cp.get_send_end());
  } else if (n == 1) {
    FAIL_IF_NOT(cp.get_num_recv_procs() == 2);
    FAIL_IF_NOT(cp.get_num_send_procs() == 2);

    // receive maps
    itr = cp.get_recv_begin();

    p = itr->first;
    cells = itr->second;

    FAIL_IF_NOT(p == 2);
    FAIL_IF_NOT(cells.size() == 4);
    FAIL_IF_NOT(cells[0] == 2);
    FAIL_IF_NOT(cells[1] == 3);
    FAIL_IF_NOT(cells[2] == 5);
    FAIL_IF_NOT(cells[3] == 6);

    itr++;

    p = itr->first;
    cells = itr->second;

    FAIL_IF_NOT(p == 3);
    FAIL_IF_NOT(cells.size() == 2);
    FAIL_IF_NOT(cells[0] == 1);
    FAIL_IF_NOT(cells[1] == 4);

    itr++;
    FAIL_IF_NOT(itr == cp.get_recv_end());

    // send maps
    itr = cp.get_send_begin();

    p = itr->first;
    cells = itr->second;

    FAIL_IF_NOT(p == 2);
    FAIL_IF_NOT(cells.size() == 1);
    FAIL_IF_NOT(cells[0] == 1);

    itr++;

    p = itr->first;
    cells = itr->second;

    FAIL_IF_NOT(p == 3);
    FAIL_IF_NOT(cells.size() == 1);
    FAIL_IF_NOT(cells[0] == 1);

    itr++;
    FAIL_IF_NOT(itr == cp.get_send_end());
  } else if (n == 2) {
    FAIL_IF_NOT(cp.get_num_recv_procs() == 3);
    FAIL_IF_NOT(cp.get_num_send_procs() == 3);

    // receive maps
    itr = cp.get_recv_begin();

    p = itr->first;
    cells = itr->second;

    FAIL_IF_NOT(p == 0);
    FAIL_IF_NOT(cells.size() == 3);
    FAIL_IF_NOT(cells[0] == 3);
    FAIL_IF_NOT(cells[1] == 5);
    FAIL_IF_NOT(cells[2] == 8);

    itr++;

    p = itr->first;
    cells = itr->second;

    FAIL_IF_NOT(p == 1);
    FAIL_IF_NOT(cells.size() == 1);
    FAIL_IF_NOT(cells[0] == 2);

    itr++;

    p = itr->first;
    cells = itr->second;

    FAIL_IF_NOT(p == 3);
    FAIL_IF_NOT(cells.size() == 4);
    FAIL_IF_NOT(cells[0] == 1);
    FAIL_IF_NOT(cells[1] == 4);
    FAIL_IF_NOT(cells[2] == 6);
    FAIL_IF_NOT(cells[3] == 7);

    itr++;
    FAIL_IF_NOT(itr == cp.get_recv_end());

    // send maps
    itr = cp.get_send_begin();

    p = itr->first;
    cells = itr->second;

    FAIL_IF_NOT(p == 0);
    FAIL_IF_NOT(cells.size() == 2);
    FAIL_IF_NOT(cells[0] == 2);
    FAIL_IF_NOT(cells[1] == 2);

    itr++;

    p = itr->first;
    cells = itr->second;

    FAIL_IF_NOT(p == 1);
    FAIL_IF_NOT(cells.size() == 4);
    FAIL_IF_NOT(cells[0] == 1);
    FAIL_IF_NOT(cells[1] == 2);
    FAIL_IF_NOT(cells[2] == 2);
    FAIL_IF_NOT(cells[3] == 1);

    itr++;

    p = itr->first;
    cells = itr->second;

    FAIL_IF_NOT(p == 3);
    FAIL_IF_NOT(cells.size() == 1);
    FAIL_IF_NOT(cells[0] == 1);

    itr++;
    FAIL_IF_NOT(itr == cp.get_send_end());
  } else if (n == 3) {
    FAIL_IF_NOT(cp.get_num_recv_procs() == 3);
    FAIL_IF_NOT(cp.get_num_send_procs() == 3);

    // receive maps
    itr = cp.get_recv_begin();

    p = itr->first;
    cells = itr->second;

    FAIL_IF_NOT(p == 0);
    FAIL_IF_NOT(cells.size() == 2);
    FAIL_IF_NOT(cells[0] == 2);
    FAIL_IF_NOT(cells[1] == 3);

    itr++;

    p = itr->first;
    cells = itr->second;

    FAIL_IF_NOT(p == 1);
    FAIL_IF_NOT(cells.size() == 1);
    FAIL_IF_NOT(cells[0] == 4);

    itr++;

    p = itr->first;
    cells = itr->second;

    FAIL_IF_NOT(p == 2);
    FAIL_IF_NOT(cells.size() == 1);
    FAIL_IF_NOT(cells[0] == 1);

    itr++;
    FAIL_IF_NOT(itr == cp.get_recv_end());

    // send maps
    itr = cp.get_send_begin();

    p = itr->first;
    cells = itr->second;

    FAIL_IF_NOT(p == 0);
    FAIL_IF_NOT(cells.size() == 2);
    FAIL_IF_NOT(cells[0] == 2);
    FAIL_IF_NOT(cells[1] == 1);

    itr++;

    p = itr->first;
    cells = itr->second;

    FAIL_IF_NOT(p == 1);
    FAIL_IF_NOT(cells.size() == 2);
    FAIL_IF_NOT(cells[0] == 1);
    FAIL_IF_NOT(cells[1] == 3);

    itr++;

    p = itr->first;
    cells = itr->second;

    FAIL_IF_NOT(p == 2);
    FAIL_IF_NOT(cells.size() == 4);
    FAIL_IF_NOT(cells[0] == 3);
    FAIL_IF_NOT(cells[1] == 1);
    FAIL_IF_NOT(cells[2] == 1);
    FAIL_IF_NOT(cells[3] == 1);

    itr++;
    FAIL_IF_NOT(itr == cp.get_send_end());
  }

  if (ut.numFails == 0) {
    ostringstream m;
    m << "DD Comm_Patterns 4 proc test with multiple "
      << "boundary cell indices ok on " << n;
    PASSMSG(m.str());
  }
}

//------------------------------------------------------------------------------------------------//
// Comm_Patterns DD Test on 9 cell mesh in autodoc/DD_Mesh.fig.  Uses the build_OS_DD_Topology()
// function from DD_Mesh.hh.

void DD_Comm_Patterns_4(rtt_dsxx::UnitTest &ut) {
  if (rtt_c4::nodes() != 4)
    return;

  // build a DD topology for a nine cell mesh
  std::shared_ptr<Topology> topology = rtt_mesh_test::build_DD_Topology();
  FAIL_IF_NOT(topology->get_parallel_scheme() == "DD");

  // make a cp object
  Comm_Patterns cp;

  // nothing should be set yet
  if (cp)
    ITFAILS;

  int n = rtt_c4::node();
  Comm_Patterns::const_iterator itr;
  int p = 0;
  vector<int> cells;

  // do three permutations of calculating the comm patterns
  for (int run = 1; run <= 3; run++) {
    // calculate comm_patterns
    cp.calc_patterns(topology);
    if (!cp)
      ITFAILS;

    // check comm_patterns
    if (n == 0) {
      FAIL_IF_NOT(cp.get_num_recv_procs() == 2);
      FAIL_IF_NOT(cp.get_num_send_procs() == 2);

      // receive maps
      itr = cp.get_recv_begin();

      p = itr->first;
      cells = itr->second;

      FAIL_IF_NOT(p == 1);
      FAIL_IF_NOT(cells.size() == 2);
      FAIL_IF_NOT(cells[0] == 1);
      FAIL_IF_NOT(cells[1] == 2);

      itr++;
      p = itr->first;
      cells = itr->second;

      FAIL_IF_NOT(p == 2);
      FAIL_IF_NOT(cells.size() == 1);
      FAIL_IF_NOT(cells[0] == 3);

      itr++;
      FAIL_IF_NOT(itr == cp.get_recv_end());

      // send maps
      itr = cp.get_send_begin();

      p = itr->first;
      cells = itr->second;

      FAIL_IF_NOT(p == 1);
      FAIL_IF_NOT(cells.size() == 2);
      FAIL_IF_NOT(cells[0] == 1);
      FAIL_IF_NOT(cells[1] == 2);

      itr++;
      p = itr->first;
      cells = itr->second;

      FAIL_IF_NOT(p == 2);
      FAIL_IF_NOT(cells.size() == 1);
      FAIL_IF_NOT(cells[0] == 2);

      itr++;
      FAIL_IF_NOT(itr == cp.get_send_end());
    } else if (n == 1) {
      FAIL_IF_NOT(cp.get_num_recv_procs() == 3);
      FAIL_IF_NOT(cp.get_num_send_procs() == 3);

      // receive maps
      itr = cp.get_recv_begin();

      p = itr->first;
      cells = itr->second;

      FAIL_IF_NOT(p == 0);
      FAIL_IF_NOT(cells.size() == 2);
      FAIL_IF_NOT(cells[0] == 1);
      FAIL_IF_NOT(cells[1] == 2);

      itr++;
      p = itr->first;
      cells = itr->second;

      FAIL_IF_NOT(p == 2);
      FAIL_IF_NOT(cells.size() == 2);
      FAIL_IF_NOT(cells[0] == 3);
      FAIL_IF_NOT(cells[1] == 4);

      itr++;
      p = itr->first;
      cells = itr->second;

      FAIL_IF_NOT(p == 3);
      FAIL_IF_NOT(cells.size() == 1);
      FAIL_IF_NOT(cells[0] == 5);

      itr++;
      FAIL_IF_NOT(itr == cp.get_recv_end());

      // send maps
      itr = cp.get_send_begin();

      p = itr->first;
      cells = itr->second;

      FAIL_IF_NOT(p == 0);
      FAIL_IF_NOT(cells.size() == 2);
      FAIL_IF_NOT(cells[0] == 1);
      FAIL_IF_NOT(cells[1] == 2);

      itr++;
      p = itr->first;
      cells = itr->second;

      FAIL_IF_NOT(p == 2);
      FAIL_IF_NOT(cells.size() == 2);
      FAIL_IF_NOT(cells[0] == 1);
      FAIL_IF_NOT(cells[1] == 2);

      itr++;
      p = itr->first;
      cells = itr->second;

      FAIL_IF_NOT(p == 3);
      FAIL_IF_NOT(cells.size() == 1);
      FAIL_IF_NOT(cells[0] == 2);

      itr++;
      FAIL_IF_NOT(itr == cp.get_send_end());
    } else if (n == 2) {
      FAIL_IF_NOT(cp.get_num_recv_procs() == 3);
      FAIL_IF_NOT(cp.get_num_send_procs() == 3);

      // receive maps
      itr = cp.get_recv_begin();

      p = itr->first;
      cells = itr->second;

      FAIL_IF_NOT(p == 0);
      FAIL_IF_NOT(cells.size() == 1);
      FAIL_IF_NOT(cells[0] == 1);

      itr++;
      p = itr->first;
      cells = itr->second;

      FAIL_IF_NOT(p == 1);
      FAIL_IF_NOT(cells.size() == 2);
      FAIL_IF_NOT(cells[0] == 2);
      FAIL_IF_NOT(cells[1] == 3);

      itr++;
      p = itr->first;
      cells = itr->second;

      FAIL_IF_NOT(p == 3);
      FAIL_IF_NOT(cells.size() == 2);
      FAIL_IF_NOT(cells[0] == 4);
      FAIL_IF_NOT(cells[1] == 5);

      itr++;
      FAIL_IF_NOT(itr == cp.get_recv_end());

      // send maps
      itr = cp.get_send_begin();

      p = itr->first;
      cells = itr->second;

      FAIL_IF_NOT(p == 0);
      FAIL_IF_NOT(cells.size() == 1);
      FAIL_IF_NOT(cells[0] == 1);

      itr++;
      p = itr->first;
      cells = itr->second;

      FAIL_IF_NOT(p == 1);
      FAIL_IF_NOT(cells.size() == 2);
      FAIL_IF_NOT(cells[0] == 1);
      FAIL_IF_NOT(cells[1] == 2);

      itr++;
      p = itr->first;
      cells = itr->second;

      FAIL_IF_NOT(p == 3);
      FAIL_IF_NOT(cells.size() == 2);
      FAIL_IF_NOT(cells[0] == 1);
      FAIL_IF_NOT(cells[1] == 2);

      itr++;
      FAIL_IF_NOT(itr == cp.get_send_end());
    } else if (n == 3) {
      FAIL_IF_NOT(cp.get_num_recv_procs() == 2);
      FAIL_IF_NOT(cp.get_num_send_procs() == 2);

      // receive maps
      itr = cp.get_recv_begin();

      p = itr->first;
      cells = itr->second;

      FAIL_IF_NOT(p == 1);
      FAIL_IF_NOT(cells.size() == 1);
      FAIL_IF_NOT(cells[0] == 1);

      itr++;
      p = itr->first;
      cells = itr->second;

      FAIL_IF_NOT(p == 2);
      FAIL_IF_NOT(cells.size() == 2);
      FAIL_IF_NOT(cells[0] == 2);
      FAIL_IF_NOT(cells[1] == 3);

      itr++;
      FAIL_IF_NOT(itr == cp.get_recv_end());

      // send maps
      itr = cp.get_send_begin();

      p = itr->first;
      cells = itr->second;

      FAIL_IF_NOT(p == 1);
      FAIL_IF_NOT(cells.size() == 1);
      FAIL_IF_NOT(cells[0] == 1);

      itr++;
      p = itr->first;
      cells = itr->second;

      FAIL_IF_NOT(p == 2);
      FAIL_IF_NOT(cells.size() == 2);
      FAIL_IF_NOT(cells[0] == 2);
      FAIL_IF_NOT(cells[1] == 3);

      itr++;
      if (itr != cp.get_send_end())
        ITFAILS;
    } else {
      ITFAILS;
    }
  }

  if (ut.numFails == 0) {
    ostringstream m;
    m << "DD Comm_Patterns 4 proc test ok on " << n;
    PASSMSG(m.str());
  }
}

//------------------------------------------------------------------------------------------------//
// Comm_Patterns full replication test.

void rep_Comm_Patterns(rtt_dsxx::UnitTest &ut) {
  // make a full replication topology for nine cell mesh
  std::shared_ptr<Topology> topology(new Rep_Topology(9));

  // make a comm_pattern object
  Comm_Patterns cp;

  if (cp)
    ITFAILS;

  // now build comm patterns.  This should be a null operation because we are in a full replication
  // topology
  cp.calc_patterns(topology);

  if (cp)
    ITFAILS;

  if (ut.numFails == 0) {
    ostringstream m;
    m << "Replication Comm_Patterns test ok on " << rtt_c4::node();
    PASSMSG(m.str());
  }
}

//------------------------------------------------------------------------------------------------//

int main(int argc, char *argv[]) {
  rtt_c4::ParallelUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    rep_Comm_Patterns(ut);
    DD_Comm_Patterns_2(ut);
    DD_Comm_Patterns_4(ut);
    DD_Comm_Patterns_4_multiple_indices(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstComm_Patterns.cc
//------------------------------------------------------------------------------------------------//
