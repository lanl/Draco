//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   mesh/test/tstParallel_Data_Op.cc
 * \author Thomas M. Evans
 * \date   Tue Mar 15 13:15:10 2005
 * \brief  Parallel_Data_Op test.
 * \note   Copyright (C) 2023 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "Test_Topology.hh"
#include "c4/Compare.hh"
#include "c4/ParallelUnitTest.hh"
#include "ds++/Release.hh"
#include "ds++/Soft_Equivalence.hh"
#include "mesh/DD_Topology.hh"
#include "mesh/Parallel_Data_Operator.hh"
#include "mesh/Rep_Topology.hh"
#include <deque>
#include <sstream>

using rtt_mesh::Comm_Patterns;
using rtt_mesh::DD_Topology;
using rtt_mesh::Parallel_Data_Operator;
using rtt_mesh::Rep_Topology;
using rtt_mesh::Topology;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

double const local_eps = std::numeric_limits<double>::epsilon();
// ------------------------------------------------------------------------------------------------ //
// TEST LOCAL-GLOBAL MAPPING
//------------------------------------------------------------------------------------------------//

void test_mapping_replication(rtt_dsxx::UnitTest &ut) {
  size_t ncells = 6;
  // build a full replication topology
  std::shared_ptr<Topology> topology(new Rep_Topology(ncells));

  // build a Parallel_Data_Operator
  Parallel_Data_Operator pdop(topology);

  // make a MT::CCSF field
  std::vector<int> local(ncells);
  std::vector<int> global(ncells);
  for (size_t cell = 0; cell < ncells; cell++)
    local[cell] = 10;

  // Do a local_global mapping with data-replicated data
  pdop.local_to_global(local, global, Parallel_Data_Operator::Data_Replicated());

  for (size_t cell = 0; cell < ncells; cell++)
    if (global[cell] != local[cell])
      ITFAILS;

  // Do a local_global mapping with data-decomposed data
  pdop.local_to_global(local, global, Parallel_Data_Operator::Data_Decomposed());

  for (size_t cell = 0; cell < ncells; cell++)
    if (global[cell] != local[cell] * rtt_c4::nodes())
      ITFAILS;

  // check local data
  for (size_t cell = 1; cell < ncells; cell++)
    if (local[cell] != 10)
      ITFAILS;

  if (ut.numFails == 0) {
    std::ostringstream m;
    m << "test_mapping_replication finished on " << rtt_c4::node();
    PASSMSG(m.str());
  }
}

//------------------------------------------------------------------------------------------------//
// DD Mapping test uses a simple topology: 2-D 6 cell mesh, 1-3 on processor 0; 4-6 on processor 1

void test_mapping_DD(rtt_dsxx::UnitTest &ut) {
  // only perform this test on two processors
  if (rtt_c4::nodes() != 2)
    return;

  // num_global = 6
  DD_Topology::Vec_Int bc(3);
  DD_Topology::Vec_Int bp(3);

  DD_Topology::Vec_Int cpp(2, 3);
  size_t global_ncells = 6;

  if (rtt_c4::node() == 0) {
    bc[0] = 1;
    bc[1] = 2;
    bc[2] = 3;

    bp[0] = 1;
    bp[1] = 1;
    bp[2] = 1;
  } else if (rtt_c4::node() == 1) {
    bc[0] = 1;
    bc[1] = 2;
    bc[2] = 3;

    bp[0] = 0;
    bp[1] = 0;
    bp[2] = 0;
  }

  // build topology
  std::shared_ptr<Topology> topology(new DD_Topology(cpp, bc, bp));

  // build a Parallel_Data_Operator
  Parallel_Data_Operator pdop(topology);

  // make fields
  std::vector<int> local(3, rtt_c4::node() + 1);
  std::vector<int> global_vect(6);

  pdop.local_to_global(local, global_vect, Parallel_Data_Operator::Data_Distributed());

  // check value of global data
  if (global_vect[0] != 1)
    ITFAILS;
  if (global_vect[1] != 1)
    ITFAILS;
  if (global_vect[2] != 1)
    ITFAILS;
  if (global_vect[3] != 2)
    ITFAILS;
  if (global_vect[4] != 2)
    ITFAILS;
  if (global_vect[5] != 2)
    ITFAILS;

  // make sure local data unchanged
  for (int loc : local)
    FAIL_IF_NOT(loc == rtt_c4::node() + 1);

  // now check mapping with different field types (vector and ccsf)
  std::vector<int> global(global_ncells);

  pdop.local_to_global(local, global, Parallel_Data_Operator::Data_Distributed());

  // check value of global data
  if (global[0] != 1)
    ITFAILS;
  if (global[1] != 1)
    ITFAILS;
  if (global[2] != 1)
    ITFAILS;
  if (global[3] != 2)
    ITFAILS;
  if (global[4] != 2)
    ITFAILS;
  if (global[5] != 2)
    ITFAILS;

  // make sure local data unchanged
  for (int loc : local)
    FAIL_IF_NOT(loc == rtt_c4::node() + 1);

  if (ut.numFails == 0) {
    std::ostringstream m;
    m << "test_mapping_DD finished on " << rtt_c4::node();
    PASSMSG(m.str());
  }
}

//------------------------------------------------------------------------------------------------//
// TEST GATHERS
//------------------------------------------------------------------------------------------------//

void gather(rtt_dsxx::UnitTest &ut) {
  if (rtt_c4::nodes() != 4)
    return;

  // make topology of 9 cell mesh
  std::shared_ptr<Topology> topology = rtt_mesh_test::build_DD_Topology();

  // make a comm patterns
  std::shared_ptr<Comm_Patterns> pattern(new Comm_Patterns());
  pattern->calc_patterns(topology);

  // make the parallel data operator class
  Parallel_Data_Operator pdop(topology);

  // make local field data on each processor
  std::vector<double> local_field(topology->num_local_cells());
  for (double &lf : local_field)
    lf = rtt_c4::node() + 10.5;

  // Get boundary cell data in a plain array:
  auto *bc_data = new double[topology->num_ghost_cells()];

  pdop.gather_ghost_cell_data(pattern, local_field, bc_data);

  // Check boundary cell data
  if (rtt_c4::node() == 0) {
    if (!rtt_dsxx::soft_equiv(bc_data[0], 11.5, local_eps))
      ITFAILS;
    if (!rtt_dsxx::soft_equiv(bc_data[1], 11.5, local_eps))
      ITFAILS;
    if (!rtt_dsxx::soft_equiv(bc_data[2], 12.5, local_eps))
      ITFAILS;
  } else if (rtt_c4::node() == 1) {
    if (!rtt_dsxx::soft_equiv(bc_data[0], 10.5, local_eps))
      ITFAILS;
    if (!rtt_dsxx::soft_equiv(bc_data[1], 10.5, local_eps))
      ITFAILS;
    if (!rtt_dsxx::soft_equiv(bc_data[2], 12.5, local_eps))
      ITFAILS;
    if (!rtt_dsxx::soft_equiv(bc_data[3], 12.5, local_eps))
      ITFAILS;
    if (!rtt_dsxx::soft_equiv(bc_data[4], 13.5, local_eps))
      ITFAILS;
  } else if (rtt_c4::node() == 2) {
    if (!rtt_dsxx::soft_equiv(bc_data[0], 10.5, local_eps))
      ITFAILS;
    if (!rtt_dsxx::soft_equiv(bc_data[1], 11.5, local_eps))
      ITFAILS;
    if (!rtt_dsxx::soft_equiv(bc_data[2], 11.5, local_eps))
      ITFAILS;
    if (!rtt_dsxx::soft_equiv(bc_data[3], 13.5, local_eps))
      ITFAILS;
    if (!rtt_dsxx::soft_equiv(bc_data[4], 13.5, local_eps))
      ITFAILS;
  } else if (rtt_c4::node() == 3) {
    if (!rtt_dsxx::soft_equiv(bc_data[0], 11.5, local_eps))
      ITFAILS;
    if (!rtt_dsxx::soft_equiv(bc_data[1], 12.5, local_eps))
      ITFAILS;
    if (!rtt_dsxx::soft_equiv(bc_data[2], 12.5, local_eps))
      ITFAILS;
  } else {
    ITFAILS;
  }

  // Get the boundary cell data in a vector
  std::vector<double> bc_data_vec;
  pdop.gather_ghost_cell_data(pattern, local_field, bc_data_vec);

  // Check the sizes of the vectors
  const size_t boundary_cells = topology->num_ghost_cells();
  if (bc_data_vec.size() != boundary_cells)
    ITFAILS;

  // Compare the results to the plain data.
  if (!rtt_dsxx::soft_equiv(bc_data_vec.begin(), bc_data_vec.end(), bc_data,
                            bc_data + boundary_cells))
    ITFAILS;

  if (ut.numFails == 0) {
    std::ostringstream m;
    m << "test_gather finished on " << rtt_c4::node();
    PASSMSG(m.str());
  }

  delete[] bc_data;
}

//------------------------------------------------------------------------------------------------//
// TEST GATHER VECTOR DATA
//------------------------------------------------------------------------------------------------//

void gather_vector_data(rtt_dsxx::UnitTest &ut) {
  if (rtt_c4::nodes() != 4)
    return;

  // make topology of 9 cell mesh
  std::shared_ptr<Topology> topology = rtt_mesh_test::build_DD_Topology();

  // make a comm patterns
  std::shared_ptr<Comm_Patterns> pattern(new Comm_Patterns());
  pattern->calc_patterns(topology);

  // make the parallel data operator class
  Parallel_Data_Operator pdop(topology);

  const size_t Ng = 10;

  // make local field data on each processor
  std::vector<std::vector<bool>> local_field(topology->num_local_cells(),
                                             std::vector<bool>(Ng, false));
  std::vector<std::vector<bool>> bc_data(topology->num_ghost_cells(), std::vector<bool>(Ng, false));
  for (auto &loc : local_field)
    for (size_t g = 0; g < Ng; g++)
      if (g < static_cast<size_t>(rtt_c4::node()))
        loc[g] = true;

  pdop.gather_ghost_cell_data(pattern, local_field, bc_data.begin(), Ng);
  if (bc_data.size() != topology->num_ghost_cells())
    ITFAILS;
  for (size_t cell = 0; cell < topology->num_ghost_cells(); cell++)
    if (bc_data[cell].size() != Ng)
      ITFAILS;

  // Check boundary cell data
  if (rtt_c4::node() == 0) {
    for (size_t g = 0; g < Ng; g++) {
      if (g < 1 && !bc_data[0][g])
        ITFAILS;
      else if (g >= 1 && bc_data[0][g])
        ITFAILS;
      if (g < 1 && !bc_data[1][g])
        ITFAILS;
      else if (g >= 1 && bc_data[1][g])
        ITFAILS;
      if (g < 2 && !bc_data[2][g])
        ITFAILS;
      else if (g >= 2 && bc_data[2][g])
        ITFAILS;
    }
  } else if (rtt_c4::node() == 1) {
    for (size_t g = 0; g < Ng; g++) {
      if (bc_data[0][g])
        ITFAILS;
      if (bc_data[1][g])
        ITFAILS;
      if (g < 2 && !bc_data[2][g])
        ITFAILS;
      else if (g >= 2 && bc_data[2][g])
        ITFAILS;
      if (g < 2 && !bc_data[3][g])
        ITFAILS;
      else if (g >= 2 && bc_data[3][g])
        ITFAILS;
      if (g < 3 && !bc_data[4][g])
        ITFAILS;
      else if (g >= 3 && bc_data[4][g])
        ITFAILS;
    }
  } else if (rtt_c4::node() == 2) {
    for (size_t g = 0; g < Ng; g++) {
      if (bc_data[0][g])
        ITFAILS;
      if (g < 1 && !bc_data[1][g])
        ITFAILS;
      else if (g >= 1 && bc_data[1][g])
        ITFAILS;
      if (g < 1 && !bc_data[2][g])
        ITFAILS;
      else if (g >= 1 && bc_data[2][g])
        ITFAILS;
      if (g < 3 && !bc_data[3][g])
        ITFAILS;
      else if (g >= 3 && bc_data[3][g])
        ITFAILS;
      if (g < 3 && !bc_data[4][g])
        ITFAILS;
      else if (g >= 3 && bc_data[4][g])
        ITFAILS;
    }
  } else if (rtt_c4::node() == 3) {
    for (size_t g = 0; g < Ng; g++) {
      if (g < 1 && !bc_data[0][g])
        ITFAILS;
      else if (g >= 1 && bc_data[0][g])
        ITFAILS;
      if (g < 2 && !bc_data[1][g])
        ITFAILS;
      else if (g >= 2 && bc_data[1][g])
        ITFAILS;
      if (g < 2 && !bc_data[2][g])
        ITFAILS;
      else if (g >= 2 && bc_data[2][g])
        ITFAILS;
    }
  } else {
    ITFAILS;
  }

  // Get the boundary cell data in a vector
  std::vector<std::vector<bool>> bc_data_vec;
  pdop.gather_ghost_cell_data(pattern, local_field, bc_data_vec, Ng);

  // Check the sizes of the vectors
  const size_t boundary_cells = topology->num_ghost_cells();
  if (bc_data_vec.size() != boundary_cells)
    ITFAILS;

  // Compare the results to the plain data.
  for (size_t cell = 0; cell < boundary_cells; cell++) {
    if (!std::equal(bc_data_vec[cell].begin(), bc_data_vec[cell].end(), bc_data[cell].begin()))
      ITFAILS;
  }

  if (ut.numFails == 0) {
    std::ostringstream m;
    m << "test_gather_vector_data finished on " << rtt_c4::node();
    PASSMSG(m.str());
  }
}

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  rtt_c4::ParallelUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    // test local-global cell mapping
    test_mapping_replication(ut);
    test_mapping_DD(ut);

    // test gathers
    gather(ut);
    gather_vector_data(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstParallel_Data_Op.cc
//------------------------------------------------------------------------------------------------//
