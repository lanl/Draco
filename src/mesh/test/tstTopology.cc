//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   mesh/test/tstTopology.cc
 * \author H. Park
 * \brief  Topology test.
 * \note   Copyright (C) 2023-2026 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "dsxx/Assert.hh"
#include "dsxx/Release.hh"
#include "dsxx/UnitTest.hh"
#include "c4/C4_Functions.hh"
#include "c4/ParallelUnitTest.hh"
#include "mesh/Comm_Patterns.hh"
#include "mesh/DD_Topology.hh"
#include "mesh/Parallel_Data_Operator.hh"
#include "mesh/Rep_Topology.hh"
#include "mesh/Topology.hh"
#include <cstddef>
#include <iostream>
#include <memory>
#include <sstream>
#include <vector>

namespace {

using rtt_mesh::DD_Topology;
using rtt_mesh::Rep_Topology;
using rtt_mesh::Topology;

//------------------------------------------------------------------------------------------------//
// HELPERS
//------------------------------------------------------------------------------------------------//

std::shared_ptr<DD_Topology> build_DD(int n) {
  // assume a 4x4 ordered mesh with cells 1-7 on proc 0 and 8-16 on proc 1

  // number of global cells; num_global = 16;
  DD_Topology::Vec_Int bc;
  DD_Topology::Vec_Int bp;

  DD_Topology::Vec_Int cpp(2);
  cpp[0] = 7;
  cpp[1] = 9;

  if (n == 0) {
    // 5 boundary cells
    bc.resize(4);
    bp.resize(4, 1);

    bc[0] = 2;
    bc[1] = 3;
    bc[2] = 4;
    bc[3] = 1;
  } else {
    // 5 boundary cells
    bc.resize(4);
    bp.resize(4, 0);

    bc[0] = 5;
    bc[1] = 6;
    bc[2] = 7;
    bc[3] = 4;
  }

  std::shared_ptr<DD_Topology> top(new DD_Topology(cpp, bc, bp));
  return top;
}

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

void test_Replication(rtt_dsxx::UnitTest &ut) {
  // assume a 20 cell mesh
  std::shared_ptr<Topology> const top(new Rep_Topology(20));
  std::shared_ptr<Rep_Topology> utop;

  FAIL_IF(top->get_parallel_scheme() != "replication");
  FAIL_IF(top->num_cells() != 20);
  FAIL_IF(top->num_local_cells() != 20);

  for (int i = 1; i <= top->num_cells(); ++i) {
    FAIL_IF(top->global_cell(i) != i);
    FAIL_IF(top->local_cell(i) != i);
  }

  {
    std::vector<char> const packed = top->pack();
    utop = std::make_shared<Rep_Topology>(packed);
    utop->print(std::cout);
  }

  const Topology &T = *utop;

  FAIL_IF(T.get_parallel_scheme() != "replication");
  FAIL_IF(T.num_cells() != 20);
  FAIL_IF(T.num_local_cells() != 20);

  for (int i = 1; i <= T.num_cells(); ++i) {
    FAIL_IF(T.global_cell(i) != i);
    FAIL_IF(T.local_cell(i) != i);
  }

  FAIL_IF(!T.valid_local_cell(1));
  FAIL_IF(!T.valid_local_cell(20));
  FAIL_IF(T.valid_local_cell(0));
  FAIL_IF(T.valid_local_cell(21));

  FAIL_IF(!T.valid_global_cell(1));
  FAIL_IF(!T.valid_global_cell(20));
  FAIL_IF(T.valid_global_cell(0));
  FAIL_IF(T.valid_global_cell(21));

  if (ut.numFails == 0)
    PASSMSG("Topology replication tests pass.");
  return;
}

//------------------------------------------------------------------------------------------------//
void test_DD(rtt_dsxx::UnitTest &ut) {
  // only perform this test on two processors
  Insist(rtt_c4::nodes() == 2, "This test only works for exactly 2 MPI ranks.");

  std::shared_ptr<Topology> top;

  if (rtt_c4::node() == 0) {
    // send off processor 1 topology
    top = build_DD(1);
    {
      std::vector<char> p = top->pack();
      int const s = static_cast<int>(p.size());
      rtt_c4::send(&s, 1, 1, 100);
      rtt_c4::send(&p[0], s, 1, 101);
    }

    // make the processor 0 topology
    top = build_DD(0);

    FAIL_IF(top->num_local_cells() != 7);
    FAIL_IF(top->global_cell(1) != 1);
    FAIL_IF(top->global_cell(2) != 2);
    FAIL_IF(top->global_cell(3) != 3);
    FAIL_IF(top->global_cell(4) != 4);
    FAIL_IF(top->global_cell(5) != 5);
    FAIL_IF(top->global_cell(6) != 6);
    FAIL_IF(top->global_cell(7) != 7);

    FAIL_IF(top->local_cell(1) != 1);
    FAIL_IF(top->local_cell(2) != 2);
    FAIL_IF(top->local_cell(3) != 3);
    FAIL_IF(top->local_cell(4) != 4);
    FAIL_IF(top->local_cell(5) != 5);
    FAIL_IF(top->local_cell(6) != 6);
    FAIL_IF(top->local_cell(7) != 7);

    for (int i = 8; i <= 16; ++i)
      FAIL_IF(top->local_cell(i) != 0);

    FAIL_IF(top->remote_local_cell(1) != 2);
    FAIL_IF(top->remote_local_cell(2) != 3);
    FAIL_IF(top->remote_local_cell(3) != 4);
    FAIL_IF(top->remote_local_cell(4) != 1);

    for (int i = 1; i <= 4; ++i)
      FAIL_IF(top->remote_proc(i) != 1);

    FAIL_IF(top->valid_local_cell(0));
    FAIL_IF(!top->valid_local_cell(5));
    FAIL_IF(top->valid_local_cell(8));
  } else {
    Check(rtt_c4::node() == 1);

    // receive topology from 0
    {
      int s = 0;
      rtt_c4::receive(&s, 1, 0, 100);
      std::vector<char> p(s);
      rtt_c4::receive(&p[0], s, 0, 101);
      top = std::make_shared<DD_Topology>(p);
    }

    FAIL_IF(top->num_local_cells() != 9);
    FAIL_IF(top->global_cell(1) != 8);
    FAIL_IF(top->global_cell(2) != 9);
    FAIL_IF(top->global_cell(3) != 10);
    FAIL_IF(top->global_cell(4) != 11);
    FAIL_IF(top->global_cell(5) != 12);
    FAIL_IF(top->global_cell(6) != 13);
    FAIL_IF(top->global_cell(7) != 14);
    FAIL_IF(top->global_cell(8) != 15);
    FAIL_IF(top->global_cell(9) != 16);

    FAIL_IF(top->local_cell(8) != 1);
    FAIL_IF(top->local_cell(9) != 2);
    FAIL_IF(top->local_cell(10) != 3);
    FAIL_IF(top->local_cell(11) != 4);
    FAIL_IF(top->local_cell(12) != 5);
    FAIL_IF(top->local_cell(13) != 6);
    FAIL_IF(top->local_cell(14) != 7);
    FAIL_IF(top->local_cell(15) != 8);
    FAIL_IF(top->local_cell(16) != 9);

    for (int i = 1; i <= 7; ++i)
      FAIL_IF(top->local_cell(i) != 0);

    FAIL_IF(top->remote_local_cell(1) != 5);
    FAIL_IF(top->remote_local_cell(2) != 6);
    FAIL_IF(top->remote_local_cell(3) != 7);
    FAIL_IF(top->remote_local_cell(4) != 4);

    for (int i = 1; i <= 4; ++i)
      FAIL_IF(top->remote_proc(i) != 0);

    FAIL_IF(top->valid_local_cell(0));
    FAIL_IF(!top->valid_local_cell(5));
    FAIL_IF(top->valid_local_cell(10));
  }
  FAIL_IF(top->get_global_offset(0) != 0);
  FAIL_IF(top->get_global_offset(1) != 7);

  FAIL_IF(top->num_cells() != 16);
  FAIL_IF(top->get_parallel_scheme() != "DD");

  FAIL_IF(top->valid_global_cell(0));
  FAIL_IF(!top->valid_global_cell(10));
  FAIL_IF(top->valid_global_cell(17));

  std::ostringstream ref_oss;
  ref_oss << "\n*** DD_TOPOLOGY DATA ***\n------------------------\n"
          << "\nParallel Topology: " << top->get_parallel_scheme() << "\n------------------------\n"
          << '\n';

  std::ostringstream test_oss;
  top->print(test_oss);
  FAIL_IF_NOT(test_oss.str() == ref_oss.str());

  std::ostringstream test2_oss;
  {
    const Topology &top_ref = *top;
    test2_oss << top_ref;
  }
  if (test2_oss.str() != test_oss.str()) {
    std::cout << test_oss.str() << '\n';
    std::cout << test2_oss.str() << '\n';
    FAILMSG("test2_oss.str() != test_oss.str().");
  }

  std::shared_ptr<rtt_mesh::Comm_Patterns> comm_pat;
  // Make the communications patterns
  comm_pat = std::make_shared<rtt_mesh::Comm_Patterns>();
  comm_pat->calc_patterns(top);

  // to identify ghost level and AMR configuration, we need to exchange level and vertex.
  std::shared_ptr<rtt_mesh::Parallel_Data_Operator> par_op;
  par_op = std::make_shared<rtt_mesh::Parallel_Data_Operator>(top);
  std::vector<int> ids(top->num_local_cells(), 0);
  for (size_t i = 0; i < ids.size(); ++i)
    ids[i] = int(i);

  par_op->fill_ghost_cell_data(comm_pat, ids);
  if (rtt_c4::node() == 0) {
    FAIL_IF(ids.size() != 11);
    FAIL_IF(ids[7] != 1);
    FAIL_IF(ids[8] != 2);
    FAIL_IF(ids[9] != 3);
    FAIL_IF(ids[10] != 0);
  } else if (rtt_c4::node() == 1) {
    FAIL_IF(ids.size() != 13);
    FAIL_IF(ids[9] != 4);
    FAIL_IF(ids[10] != 5);
    FAIL_IF(ids[11] != 6);
    FAIL_IF(ids[12] != 3);
  }

  if (ut.numFails == 0) {
    std::ostringstream m;
    m << "Topology DD tests pass on " << rtt_c4::node();
    PASSMSG(m.str());
  }
  return;
}

} // namespace

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  rtt_c4::ParallelUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    // full replication tests of General_Topology and Rep_Topology
    test_Replication(ut);

    // DD tests of General_Topology (2 PE case only)
    if (rtt_c4::nodes() == 2)
      test_DD(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstTopology.cc
//------------------------------------------------------------------------------------------------//
