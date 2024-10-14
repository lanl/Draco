//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   mesh/test/tstTopology_Builder.cc
 * \author Thomas M. Evans
 * \date   Wed May 11 09:47:27 2005
 * \brief  Topology_Builder test.
 * \note   Copyright (C) 2023-2024 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "Test_Topology.hh"
#include "dsxx/Assert.hh"
#include "dsxx/Release.hh"
#include "dsxx/UnitTest.hh"
#include "c4/C4_Functions.hh"
#include "c4/ParallelUnitTest.hh"
#include "mesh/Topology_Builder.hh"
#include <iostream>
#include <memory>
#include <sstream>

using rtt_mesh::Topology;

using TBuilder = rtt_mesh::Topology_Builder;
using Replication = rtt_mesh_test::One_Proc;
using Two_P = rtt_mesh_test::Two_Proc;
using Two_P_NGC = rtt_mesh_test::Two_Proc_NGC;
using Four_P = rtt_mesh_test::Four_Proc;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//
// tests replication of an arbitrary 6 cell topology

void replication_test(rtt_dsxx::UnitTest &ut) {
  Replication const iface;

  // make a topology builder
  TBuilder topology_builder(iface.decomposition_scheme(), iface.get_num_global_cells(),
                            iface.get_num_bnd_cells(), iface.get_bnd_cells(), iface.get_bnd_procs(),
                            iface.has_ghost_cells() ? iface.get_num_full_cells() : 0,
                            iface.get_num_cells(), iface.get_full_active());

  // get the topology
  std::shared_ptr<Topology> const topology = topology_builder.build_Topology();
  FAIL_IF(!topology);
  FAIL_IF(topology->get_parallel_scheme() != "replication");
  FAIL_IF(topology->num_cells() != 6);
  FAIL_IF(topology->num_local_cells() != 6);

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
  Two_P const iface;

  // make a topology builder
  TBuilder topology_builder(iface.decomposition_scheme(), iface.get_num_global_cells(),
                            iface.get_num_bnd_cells(), iface.get_bnd_cells(), iface.get_bnd_procs(),
                            iface.has_ghost_cells() ? iface.get_num_full_cells() : 0,
                            iface.get_num_cells(), iface.get_full_active());

  // get the topology
  std::shared_ptr<Topology> const topology = topology_builder.build_Topology();
  FAIL_IF(!topology);

  // check the topologies
  FAIL_IF(topology->get_parallel_scheme() != "DD");
  FAIL_IF(topology->num_cells() != 15);
  // check decompositions
  if (rtt_c4::node() == 0) {
    FAIL_IF(topology->num_local_cells() != 6);
    FAIL_IF(topology->num_ghost_cells() != 3);

    FAIL_IF(topology->remote_local_cell(1) != 1);
    FAIL_IF(topology->remote_local_cell(2) != 3);
    FAIL_IF(topology->remote_local_cell(3) != 5);

    FAIL_IF(topology->remote_proc(1) != 1);
    FAIL_IF(topology->remote_proc(2) != 1);
    FAIL_IF(topology->remote_proc(3) != 1);

    FAIL_IF(topology->global_cell(1) != 1);
    FAIL_IF(topology->global_cell(2) != 2);
    FAIL_IF(topology->global_cell(3) != 3);
    FAIL_IF(topology->global_cell(4) != 4);
    FAIL_IF(topology->global_cell(5) != 5);
    FAIL_IF(topology->global_cell(6) != 6);
  } else {
    FAIL_IF(topology->num_local_cells() != 9);
    FAIL_IF(topology->num_ghost_cells() != 2);

    FAIL_IF(topology->remote_local_cell(1) != 6);
    FAIL_IF(topology->remote_local_cell(2) != 5);

    FAIL_IF(topology->remote_proc(1) != 0);
    FAIL_IF(topology->remote_proc(2) != 0);

    FAIL_IF(topology->global_cell(1) != 7);
    FAIL_IF(topology->global_cell(2) != 8);
    FAIL_IF(topology->global_cell(3) != 9);
    FAIL_IF(topology->global_cell(4) != 10);
    FAIL_IF(topology->global_cell(5) != 11);
    FAIL_IF(topology->global_cell(6) != 12);
    FAIL_IF(topology->global_cell(7) != 13);
    FAIL_IF(topology->global_cell(8) != 14);
    FAIL_IF(topology->global_cell(9) != 15);
  }
  std::cout << "and now we at?\n";

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
  Two_P_NGC const iface;

  // make a topology builder
  TBuilder topology_builder(iface.decomposition_scheme(), iface.get_num_global_cells(),
                            iface.get_num_bnd_cells(), iface.get_bnd_cells(), iface.get_bnd_procs(),
                            iface.has_ghost_cells() ? iface.get_num_full_cells() : 0,
                            iface.get_num_cells(), iface.get_full_active());

  // get the topology
  std::shared_ptr<Topology> const topology = topology_builder.build_Topology();
  FAIL_IF(!topology);

  // check the topologies
  FAIL_IF(topology->get_parallel_scheme() != "DD");
  FAIL_IF(topology->num_cells() != 15);

  // check decompositions
  if (rtt_c4::node() == 0) {
    FAIL_IF(topology->num_local_cells() != 6);
    FAIL_IF(topology->num_ghost_cells() != 3);

    FAIL_IF(topology->remote_local_cell(1) != 1);
    FAIL_IF(topology->remote_local_cell(2) != 3);
    FAIL_IF(topology->remote_local_cell(3) != 5);

    FAIL_IF(topology->remote_proc(1) != 1);
    FAIL_IF(topology->remote_proc(2) != 1);
    FAIL_IF(topology->remote_proc(3) != 1);

    FAIL_IF(topology->global_cell(1) != 1);
    FAIL_IF(topology->global_cell(2) != 2);
    FAIL_IF(topology->global_cell(3) != 3);
    FAIL_IF(topology->global_cell(4) != 4);
    FAIL_IF(topology->global_cell(5) != 5);
    FAIL_IF(topology->global_cell(6) != 6);
  } else {
    FAIL_IF(topology->num_local_cells() != 9);
    FAIL_IF(topology->num_ghost_cells() != 2);

    FAIL_IF(topology->remote_local_cell(1) != 6);
    FAIL_IF(topology->remote_local_cell(2) != 5);

    FAIL_IF(topology->remote_proc(1) != 0);
    FAIL_IF(topology->remote_proc(2) != 0);

    FAIL_IF(topology->global_cell(1) != 7);
    FAIL_IF(topology->global_cell(2) != 8);
    FAIL_IF(topology->global_cell(3) != 9);
    FAIL_IF(topology->global_cell(4) != 10);
    FAIL_IF(topology->global_cell(5) != 11);
    FAIL_IF(topology->global_cell(6) != 12);
    FAIL_IF(topology->global_cell(7) != 13);
    FAIL_IF(topology->global_cell(8) != 14);
    FAIL_IF(topology->global_cell(9) != 15);
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
  Four_P const iface;

  // make a topology builder
  TBuilder topology_builder(iface.decomposition_scheme(), iface.get_num_global_cells(),
                            iface.get_num_bnd_cells(), iface.get_bnd_cells(), iface.get_bnd_procs(),
                            iface.has_ghost_cells() ? iface.get_num_full_cells() : 0,
                            iface.get_num_cells(), iface.get_full_active());

  // get the topology
  std::shared_ptr<Topology> const topology = topology_builder.build_Topology();
  FAIL_IF(!topology);

  // check the topologies
  FAIL_IF(topology->get_parallel_scheme() != "DD");
  FAIL_IF(topology->num_cells() != 18);

  // check decompositions
  if (rtt_c4::node() == 0) {
    FAIL_IF(topology->num_local_cells() != 2);
    FAIL_IF(topology->num_ghost_cells() != 8);

    FAIL_IF(topology->remote_local_cell(1) != 2);
    FAIL_IF(topology->remote_local_cell(2) != 4);
    FAIL_IF(topology->remote_local_cell(3) != 5);
    FAIL_IF(topology->remote_local_cell(4) != 1);
    FAIL_IF(topology->remote_local_cell(5) != 1);
    FAIL_IF(topology->remote_local_cell(6) != 3);
    FAIL_IF(topology->remote_local_cell(7) != 5);
    FAIL_IF(topology->remote_local_cell(8) != 7);

    FAIL_IF(topology->remote_proc(1) != 3);
    FAIL_IF(topology->remote_proc(2) != 3);
    FAIL_IF(topology->remote_proc(3) != 3);
    FAIL_IF(topology->remote_proc(4) != 1);
    FAIL_IF(topology->remote_proc(5) != 2);
    FAIL_IF(topology->remote_proc(6) != 2);
    FAIL_IF(topology->remote_proc(7) != 2);
    FAIL_IF(topology->remote_proc(8) != 2);

    FAIL_IF(topology->global_cell(1) != 1);
    FAIL_IF(topology->global_cell(2) != 2);
  }

  if (rtt_c4::node() == 1) {
    FAIL_IF(topology->num_local_cells() != 2);
    FAIL_IF(topology->num_ghost_cells() != 4);

    FAIL_IF(topology->remote_local_cell(1) != 6);
    FAIL_IF(topology->remote_local_cell(2) != 2);
    FAIL_IF(topology->remote_local_cell(3) != 7);
    FAIL_IF(topology->remote_local_cell(4) != 8);

    FAIL_IF(topology->remote_proc(1) != 3);
    FAIL_IF(topology->remote_proc(2) != 0);
    FAIL_IF(topology->remote_proc(3) != 2);
    FAIL_IF(topology->remote_proc(4) != 2);

    FAIL_IF(topology->global_cell(1) != 3);
    FAIL_IF(topology->global_cell(2) != 4);
  }

  if (rtt_c4::node() == 2) {
    FAIL_IF(topology->num_local_cells() != 8);
    FAIL_IF(topology->num_ghost_cells() != 3);

    FAIL_IF(topology->remote_local_cell(1) != 1);
    FAIL_IF(topology->remote_local_cell(2) != 2);
    FAIL_IF(topology->remote_local_cell(3) != 2);

    FAIL_IF(topology->remote_proc(1) != 0);
    FAIL_IF(topology->remote_proc(2) != 0);
    FAIL_IF(topology->remote_proc(3) != 1);

    FAIL_IF(topology->global_cell(1) != 5);
    FAIL_IF(topology->global_cell(2) != 6);
    FAIL_IF(topology->global_cell(3) != 7);
    FAIL_IF(topology->global_cell(4) != 8);
    FAIL_IF(topology->global_cell(5) != 9);
    FAIL_IF(topology->global_cell(6) != 10);
    FAIL_IF(topology->global_cell(7) != 11);
    FAIL_IF(topology->global_cell(8) != 12);
  }

  if (rtt_c4::node() == 3) {
    FAIL_IF(topology->num_local_cells() != 6);
    FAIL_IF(topology->num_ghost_cells() != 3);

    FAIL_IF(topology->remote_local_cell(1) != 1);
    FAIL_IF(topology->remote_local_cell(2) != 2);
    FAIL_IF(topology->remote_local_cell(3) != 1);

    FAIL_IF(topology->remote_proc(1) != 0);
    FAIL_IF(topology->remote_proc(2) != 0);
    FAIL_IF(topology->remote_proc(3) != 1);

    FAIL_IF(topology->global_cell(1) != 13);
    FAIL_IF(topology->global_cell(2) != 14);
    FAIL_IF(topology->global_cell(3) != 15);
    FAIL_IF(topology->global_cell(4) != 16);
    FAIL_IF(topology->global_cell(5) != 17);
    FAIL_IF(topology->global_cell(6) != 18);
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
