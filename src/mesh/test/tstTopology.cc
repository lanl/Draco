//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   mesh/test/tstTopology.cc
 * \author H. Park
 * \date   Aug 1 2023
 * \brief  Topology test.
 * \note   Copyright (C) 2023 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "c4/ParallelUnitTest.hh"
#include "ds++/Release.hh"
#include "mesh/DD_Topology.hh"
#include "mesh/Rep_Topology.hh"
#include <sstream>

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
  std::shared_ptr<Topology> top(new Rep_Topology(20));
  std::shared_ptr<Rep_Topology> utop;

  if (top->get_parallel_scheme() != "replication")
    ITFAILS;
  if (top->num_cells() != 20)
    ITFAILS;
  if (top->num_local_cells() != 20)
    ITFAILS;

  for (int i = 1; i <= top->num_cells(); ++i) {
    if (top->global_cell(i) != i)
      ITFAILS;
    if (top->local_cell(i) != i)
      ITFAILS;
  }

  {
    std::vector<char> packed = top->pack();
    utop = std::make_shared<Rep_Topology>(packed);
    utop->print(std::cout);
  }

  const Topology &T = *utop;

  if (T.get_parallel_scheme() != "replication")
    ITFAILS;
  if (T.num_cells() != 20)
    ITFAILS;
  if (T.num_local_cells() != 20)
    ITFAILS;

  for (int i = 1; i <= T.num_cells(); ++i) {
    if (T.global_cell(i) != i)
      ITFAILS;
    if (T.local_cell(i) != i)
      ITFAILS;
  }

  if (!T.valid_local_cell(1))
    ITFAILS;
  if (!T.valid_local_cell(20))
    ITFAILS;
  if (T.valid_local_cell(0))
    ITFAILS;
  if (T.valid_local_cell(21))
    ITFAILS;

  if (!T.valid_global_cell(1))
    ITFAILS;
  if (!T.valid_global_cell(20))
    ITFAILS;
  if (T.valid_global_cell(0))
    ITFAILS;
  if (T.valid_global_cell(21))
    ITFAILS;

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
      int s = static_cast<int>(p.size());
      rtt_c4::send(&s, 1, 1, 100);
      rtt_c4::send(&p[0], s, 1, 101);
    }

    // make the processor 0 topology
    top = build_DD(0);

    if (top->num_local_cells() != 7)
      ITFAILS;
    if (top->global_cell(1) != 1)
      ITFAILS;
    if (top->global_cell(2) != 2)
      ITFAILS;
    if (top->global_cell(3) != 3)
      ITFAILS;
    if (top->global_cell(4) != 4)
      ITFAILS;
    if (top->global_cell(5) != 5)
      ITFAILS;
    if (top->global_cell(6) != 6)
      ITFAILS;
    if (top->global_cell(7) != 7)
      ITFAILS;

    if (top->local_cell(1) != 1)
      ITFAILS;
    if (top->local_cell(2) != 2)
      ITFAILS;
    if (top->local_cell(3) != 3)
      ITFAILS;
    if (top->local_cell(4) != 4)
      ITFAILS;
    if (top->local_cell(5) != 5)
      ITFAILS;
    if (top->local_cell(6) != 6)
      ITFAILS;
    if (top->local_cell(7) != 7)
      ITFAILS;

    for (int i = 8; i <= 16; ++i)
      if (top->local_cell(i) != 0)
        ITFAILS;

    if (top->remote_local_cell(1) != 2)
      ITFAILS;
    if (top->remote_local_cell(2) != 3)
      ITFAILS;
    if (top->remote_local_cell(3) != 4)
      ITFAILS;
    if (top->remote_local_cell(4) != 1)
      ITFAILS;

    for (int i = 1; i <= 4; ++i)
      if (top->remote_proc(i) != 1)
        ITFAILS;

    if (top->valid_local_cell(0))
      ITFAILS;
    if (!top->valid_local_cell(5))
      ITFAILS;
    if (top->valid_local_cell(8))
      ITFAILS;
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

    if (top->num_local_cells() != 9)
      ITFAILS;
    if (top->global_cell(1) != 8)
      ITFAILS;
    if (top->global_cell(2) != 9)
      ITFAILS;
    if (top->global_cell(3) != 10)
      ITFAILS;
    if (top->global_cell(4) != 11)
      ITFAILS;
    if (top->global_cell(5) != 12)
      ITFAILS;
    if (top->global_cell(6) != 13)
      ITFAILS;
    if (top->global_cell(7) != 14)
      ITFAILS;
    if (top->global_cell(8) != 15)
      ITFAILS;
    if (top->global_cell(9) != 16)
      ITFAILS;

    if (top->local_cell(8) != 1)
      ITFAILS;
    if (top->local_cell(9) != 2)
      ITFAILS;
    if (top->local_cell(10) != 3)
      ITFAILS;
    if (top->local_cell(11) != 4)
      ITFAILS;
    if (top->local_cell(12) != 5)
      ITFAILS;
    if (top->local_cell(13) != 6)
      ITFAILS;
    if (top->local_cell(14) != 7)
      ITFAILS;
    if (top->local_cell(15) != 8)
      ITFAILS;
    if (top->local_cell(16) != 9)
      ITFAILS;

    for (int i = 1; i <= 7; ++i)
      if (top->local_cell(i) != 0)
        ITFAILS;

    if (top->remote_local_cell(1) != 5)
      ITFAILS;
    if (top->remote_local_cell(2) != 6)
      ITFAILS;
    if (top->remote_local_cell(3) != 7)
      ITFAILS;
    if (top->remote_local_cell(4) != 4)
      ITFAILS;

    for (int i = 1; i <= 4; ++i)
      if (top->remote_proc(i) != 0)
        ITFAILS;

    if (top->valid_local_cell(0))
      ITFAILS;
    if (!top->valid_local_cell(5))
      ITFAILS;
    if (top->valid_local_cell(10))
      ITFAILS;
  }
  if (top->get_global_offset(0) != 0)
    ITFAILS;
  if (top->get_global_offset(1) != 7)
    ITFAILS;

  if (top->num_cells() != 16)
    ITFAILS;
  if (top->get_parallel_scheme() != "DD")
    ITFAILS;

  if (top->valid_global_cell(0))
    ITFAILS;
  if (!top->valid_global_cell(10))
    ITFAILS;
  if (top->valid_global_cell(17))
    ITFAILS;

  std::ostringstream ref_oss;
  ref_oss << "\n*** DD_TOPOLOGY DATA ***"
          << "\n------------------------\n"
          << "\nParallel Topology: " << top->get_parallel_scheme() << "\n------------------------\n"
          << std::endl;

  std::ostringstream test_oss;
  top->print(test_oss);
  if (test_oss.str() != ref_oss.str())
    ITFAILS;

  std::ostringstream test2_oss;
  {
    const Topology &top_ref = *top;
    test2_oss << top_ref;
  }
  if (test2_oss.str() != test_oss.str()) {
    std::cout << test_oss.str() << std::endl;
    std::cout << test2_oss.str() << std::endl;
    ITFAILS;
  }

  if (ut.numFails == 0) {
    std::ostringstream m;
    m << "Topology DD tests pass on " << rtt_c4::node();
    PASSMSG(m.str());
  }
  return;
}

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
