//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   mesh/test/tstTopology_Builder.cc
 * \author Thomas M. Evans
 * \date   Wed May 11 09:47:27 2005
 * \brief  Topology_Builder test.
 * \note   Copyright (C) 2023 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

#include "Test_Topology.hh"
#include "c4/ParallelUnitTest.hh"
#include "ds++/Release.hh"
#include "mesh/Topology_Builder.hh"
#include <sstream>
#include <vector>

using rtt_mesh::Topology;

using TBuilder = rtt_mesh::Topology_Builder;
using Replication = rtt_mesh_test::One_Proc;
using Two_P = rtt_mesh_test::Two_Proc;
using Two_P_NGC = rtt_mesh_test::Two_Proc_NGC;
using Four_P = rtt_mesh_test::Four_Proc;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//
// tests replication of an arbitary 6 cell topology

void replication_test(rtt_dsxx::UnitTest &ut) {
  Replication iface;

  // make a topology builder
  TBuilder topology_builder(iface.decomposition_scheme(), iface.get_num_global_cells(),
                            iface.get_num_bnd_cells(), iface.get_bnd_cells(), iface.get_bnd_procs(),
                            iface.has_ghost_cells() ? iface.get_num_full_cells() : 0,
                            iface.get_num_cells(), iface.get_full_active());

  // get the topology
  std::shared_ptr<Topology> topology = topology_builder.build_Topology();
  if (!topology)
    ITFAILS;

  if (topology->get_parallel_scheme() != "replication")
    ITFAILS;
  if (topology->num_cells() != 6)
    ITFAILS;
  if (topology->num_local_cells() != 6)
    ITFAILS;

  if (ut.numFails == 0) {
    std::ostringstream m;
    m << "Replication test passed on " << rtt_c4::node();
    PASSMSG(m.str());
  }
}

//------------------------------------------------------------------------------------------------//
// uses the 15 cell mesh in autodoc/2_proc_DD.fig

void test_2_proc(rtt_dsxx::UnitTest &ut) {
  Require(rtt_c4::nodes() == 2);

  // make a 2 proc interface
  Two_P iface;

  // make a topology builder
  TBuilder topology_builder(iface.decomposition_scheme(), iface.get_num_global_cells(),
                            iface.get_num_bnd_cells(), iface.get_bnd_cells(), iface.get_bnd_procs(),
                            iface.has_ghost_cells() ? iface.get_num_full_cells() : 0,
                            iface.get_num_cells(), iface.get_full_active());

  // get the topology
  std::shared_ptr<Topology> topology = topology_builder.build_Topology();
  if (!topology)
    ITFAILS;

  // check the topologies
  if (topology->get_parallel_scheme() != "DD")
    ITFAILS;
  if (topology->num_cells() != 15)
    ITFAILS;
  // check decompositions
  if (rtt_c4::node() == 0) {
    if (topology->num_local_cells() != 6)
      ITFAILS;
    if (topology->num_ghost_cells() != 3)
      ITFAILS;

    if (topology->remote_local_cell(1) != 1)
      ITFAILS;
    if (topology->remote_local_cell(2) != 3)
      ITFAILS;
    if (topology->remote_local_cell(3) != 5)
      ITFAILS;

    if (topology->remote_proc(1) != 1)
      ITFAILS;
    if (topology->remote_proc(2) != 1)
      ITFAILS;
    if (topology->remote_proc(3) != 1)
      ITFAILS;

    if (topology->global_cell(1) != 1)
      ITFAILS;
    if (topology->global_cell(2) != 2)
      ITFAILS;
    if (topology->global_cell(3) != 3)
      ITFAILS;
    if (topology->global_cell(4) != 4)
      ITFAILS;
    if (topology->global_cell(5) != 5)
      ITFAILS;
    if (topology->global_cell(6) != 6)
      ITFAILS;
  } else {
    if (topology->num_local_cells() != 9)
      ITFAILS;
    if (topology->num_ghost_cells() != 2)
      ITFAILS;

    if (topology->remote_local_cell(1) != 6)
      ITFAILS;
    if (topology->remote_local_cell(2) != 5)
      ITFAILS;

    if (topology->remote_proc(1) != 0)
      ITFAILS;
    if (topology->remote_proc(2) != 0)
      ITFAILS;

    if (topology->global_cell(1) != 7)
      ITFAILS;
    if (topology->global_cell(2) != 8)
      ITFAILS;
    if (topology->global_cell(3) != 9)
      ITFAILS;
    if (topology->global_cell(4) != 10)
      ITFAILS;
    if (topology->global_cell(5) != 11)
      ITFAILS;
    if (topology->global_cell(6) != 12)
      ITFAILS;
    if (topology->global_cell(7) != 13)
      ITFAILS;
    if (topology->global_cell(8) != 14)
      ITFAILS;
    if (topology->global_cell(9) != 15)
      ITFAILS;
  }
  std::cout << "and now we at? " << std::endl;

  if (ut.numFails == 0) {
    std::ostringstream m;
    m << "2 Proc DD test passed on " << rtt_c4::node();
    PASSMSG(m.str());
  }
}

//------------------------------------------------------------------------------------------------//
// uses the 15 cell mesh in autodoc/2_proc_DD.fig without ghost cells

void test_2_proc_no_ghost_cells(rtt_dsxx::UnitTest &ut) {
  Require(rtt_c4::nodes() == 2);

  // make a 2 proc interface
  Two_P_NGC iface;

  // make a topology builder
  TBuilder topology_builder(iface.decomposition_scheme(), iface.get_num_global_cells(),
                            iface.get_num_bnd_cells(), iface.get_bnd_cells(), iface.get_bnd_procs(),
                            iface.has_ghost_cells() ? iface.get_num_full_cells() : 0,
                            iface.get_num_cells(), iface.get_full_active());

  // get the topology
  std::shared_ptr<Topology> topology = topology_builder.build_Topology();
  if (!topology)
    ITFAILS;

  // check the topologies
  if (topology->get_parallel_scheme() != "DD")
    ITFAILS;
  if (topology->num_cells() != 15)
    ITFAILS;

  // check decompositions
  if (rtt_c4::node() == 0) {
    if (topology->num_local_cells() != 6)
      ITFAILS;
    if (topology->num_ghost_cells() != 3)
      ITFAILS;

    if (topology->remote_local_cell(1) != 1)
      ITFAILS;
    if (topology->remote_local_cell(2) != 3)
      ITFAILS;
    if (topology->remote_local_cell(3) != 5)
      ITFAILS;

    if (topology->remote_proc(1) != 1)
      ITFAILS;
    if (topology->remote_proc(2) != 1)
      ITFAILS;
    if (topology->remote_proc(3) != 1)
      ITFAILS;

    if (topology->global_cell(1) != 1)
      ITFAILS;
    if (topology->global_cell(2) != 2)
      ITFAILS;
    if (topology->global_cell(3) != 3)
      ITFAILS;
    if (topology->global_cell(4) != 4)
      ITFAILS;
    if (topology->global_cell(5) != 5)
      ITFAILS;
    if (topology->global_cell(6) != 6)
      ITFAILS;
  } else {
    if (topology->num_local_cells() != 9)
      ITFAILS;
    if (topology->num_ghost_cells() != 2)
      ITFAILS;

    if (topology->remote_local_cell(1) != 6)
      ITFAILS;
    if (topology->remote_local_cell(2) != 5)
      ITFAILS;

    if (topology->remote_proc(1) != 0)
      ITFAILS;
    if (topology->remote_proc(2) != 0)
      ITFAILS;

    if (topology->global_cell(1) != 7)
      ITFAILS;
    if (topology->global_cell(2) != 8)
      ITFAILS;
    if (topology->global_cell(3) != 9)
      ITFAILS;
    if (topology->global_cell(4) != 10)
      ITFAILS;
    if (topology->global_cell(5) != 11)
      ITFAILS;
    if (topology->global_cell(6) != 12)
      ITFAILS;
    if (topology->global_cell(7) != 13)
      ITFAILS;
    if (topology->global_cell(8) != 14)
      ITFAILS;
    if (topology->global_cell(9) != 15)
      ITFAILS;
  }

  if (ut.numFails == 0) {
    std::ostringstream m;
    m << "2 Proc DD test without ghost cells passed on " << rtt_c4::node();
    PASSMSG(m.str());
  }
}

//------------------------------------------------------------------------------------------------//
// uses the 18 cell mesh in autodoc/4_proc_DD.fig

void test_4_proc(rtt_dsxx::UnitTest &ut) {
  Require(rtt_c4::nodes() == 4);

  // make a 2 proc interface
  Four_P iface;

  // make a topology builder
  TBuilder topology_builder(iface.decomposition_scheme(), iface.get_num_global_cells(),
                            iface.get_num_bnd_cells(), iface.get_bnd_cells(), iface.get_bnd_procs(),
                            iface.has_ghost_cells() ? iface.get_num_full_cells() : 0,
                            iface.get_num_cells(), iface.get_full_active());

  // get the topology
  std::shared_ptr<Topology> topology = topology_builder.build_Topology();
  if (!topology)
    ITFAILS;

  // check the topologies
  if (topology->get_parallel_scheme() != "DD")
    ITFAILS;
  if (topology->num_cells() != 18)
    ITFAILS;

  // check decompositions
  if (rtt_c4::node() == 0) {
    if (topology->num_local_cells() != 2)
      ITFAILS;
    if (topology->num_ghost_cells() != 8)
      ITFAILS;

    if (topology->remote_local_cell(1) != 2)
      ITFAILS;
    if (topology->remote_local_cell(2) != 4)
      ITFAILS;
    if (topology->remote_local_cell(3) != 5)
      ITFAILS;
    if (topology->remote_local_cell(4) != 1)
      ITFAILS;
    if (topology->remote_local_cell(5) != 1)
      ITFAILS;
    if (topology->remote_local_cell(6) != 3)
      ITFAILS;
    if (topology->remote_local_cell(7) != 5)
      ITFAILS;
    if (topology->remote_local_cell(8) != 7)
      ITFAILS;

    if (topology->remote_proc(1) != 3)
      ITFAILS;
    if (topology->remote_proc(2) != 3)
      ITFAILS;
    if (topology->remote_proc(3) != 3)
      ITFAILS;
    if (topology->remote_proc(4) != 1)
      ITFAILS;
    if (topology->remote_proc(5) != 2)
      ITFAILS;
    if (topology->remote_proc(6) != 2)
      ITFAILS;
    if (topology->remote_proc(7) != 2)
      ITFAILS;
    if (topology->remote_proc(8) != 2)
      ITFAILS;

    if (topology->global_cell(1) != 1)
      ITFAILS;
    if (topology->global_cell(2) != 2)
      ITFAILS;
  }

  if (rtt_c4::node() == 1) {
    if (topology->num_local_cells() != 2)
      ITFAILS;
    if (topology->num_ghost_cells() != 4)
      ITFAILS;

    if (topology->remote_local_cell(1) != 6)
      ITFAILS;
    if (topology->remote_local_cell(2) != 2)
      ITFAILS;
    if (topology->remote_local_cell(3) != 7)
      ITFAILS;
    if (topology->remote_local_cell(4) != 8)
      ITFAILS;

    if (topology->remote_proc(1) != 3)
      ITFAILS;
    if (topology->remote_proc(2) != 0)
      ITFAILS;
    if (topology->remote_proc(3) != 2)
      ITFAILS;
    if (topology->remote_proc(4) != 2)
      ITFAILS;

    if (topology->global_cell(1) != 3)
      ITFAILS;
    if (topology->global_cell(2) != 4)
      ITFAILS;
  }

  if (rtt_c4::node() == 2) {
    if (topology->num_local_cells() != 8)
      ITFAILS;
    if (topology->num_ghost_cells() != 3)
      ITFAILS;

    if (topology->remote_local_cell(1) != 1)
      ITFAILS;
    if (topology->remote_local_cell(2) != 2)
      ITFAILS;
    if (topology->remote_local_cell(3) != 2)
      ITFAILS;

    if (topology->remote_proc(1) != 0)
      ITFAILS;
    if (topology->remote_proc(2) != 0)
      ITFAILS;
    if (topology->remote_proc(3) != 1)
      ITFAILS;

    if (topology->global_cell(1) != 5)
      ITFAILS;
    if (topology->global_cell(2) != 6)
      ITFAILS;
    if (topology->global_cell(3) != 7)
      ITFAILS;
    if (topology->global_cell(4) != 8)
      ITFAILS;
    if (topology->global_cell(5) != 9)
      ITFAILS;
    if (topology->global_cell(6) != 10)
      ITFAILS;
    if (topology->global_cell(7) != 11)
      ITFAILS;
    if (topology->global_cell(8) != 12)
      ITFAILS;
  }

  if (rtt_c4::node() == 3) {
    if (topology->num_local_cells() != 6)
      ITFAILS;
    if (topology->num_ghost_cells() != 3)
      ITFAILS;

    if (topology->remote_local_cell(1) != 1)
      ITFAILS;
    if (topology->remote_local_cell(2) != 2)
      ITFAILS;
    if (topology->remote_local_cell(3) != 1)
      ITFAILS;

    if (topology->remote_proc(1) != 0)
      ITFAILS;
    if (topology->remote_proc(2) != 0)
      ITFAILS;
    if (topology->remote_proc(3) != 1)
      ITFAILS;

    if (topology->global_cell(1) != 13)
      ITFAILS;
    if (topology->global_cell(2) != 14)
      ITFAILS;
    if (topology->global_cell(3) != 15)
      ITFAILS;
    if (topology->global_cell(4) != 16)
      ITFAILS;
    if (topology->global_cell(5) != 17)
      ITFAILS;
    if (topology->global_cell(6) != 18)
      ITFAILS;
  }

  if (ut.numFails == 0) {
    std::ostringstream m;
    m << "4 Proc DD test passed on " << rtt_c4::node();
    PASSMSG(m.str());
  }
}

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  rtt_c4::ParallelUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    replication_test(ut);

    if (rtt_c4::nodes() == 2) {
      test_2_proc(ut);
      test_2_proc_no_ghost_cells(ut);
    }
    if (rtt_c4::nodes() == 4) {
      test_4_proc(ut);
    }
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of mesh/test/tstTopology_Builder.cc
//------------------------------------------------------------------------------------------------//
