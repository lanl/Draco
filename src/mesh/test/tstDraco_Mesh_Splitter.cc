//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   mesh/test/tstDraco_Mesh_Splitter.cc
 * \author Ryan Wollaeger <wollaeger@lanl.gov>
 * \brief  Draco_Mesh_Splitter class unit test.
 * \note   Copyright (C) 2025-2026 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "dsxx/Assert.hh"
#include "dsxx/Release.hh"
#include "dsxx/Soft_Equivalence.hh"
#include "dsxx/UnitTest.hh"
#include "c4/C4_Functions.hh"
#include "c4/ParallelUnitTest.hh"
#include "mesh/Draco_Mesh.hh"
#include "mesh/Draco_Mesh_Builder.hh"
#include "mesh/Draco_Mesh_Splitter.hh"
#include "mesh/X3D_Draco_Mesh_Reader.hh"
#include "mesh_element/Geometry.hh"
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

// cSpell:ignore bdy

namespace {

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

void split_1d_mesh_dd(rtt_c4::ParallelUnitTest &ut) {
  //----------------------------------------------------------------------------------------------//
  // Two ranks, 3 global cells, 1D mesh
  //
  // | r0 | r1 | r1 |
  // | c0 | c1 | c2 |
  // -----------------> x
  // n0   n1   n2   n3
  //----------------------------------------------------------------------------------------------//

  Insist(rtt_c4::nodes() == 2, "This test only uses 2 PE.");

  // create global mesh
  const unsigned dim = 1;
  const rtt_mesh_element::Geometry geom = rtt_mesh_element::Geometry::CARTESIAN;
  const std::vector<unsigned> num_faces_per_cell = {2, 2, 2};
  const std::vector<unsigned> cell_to_node_linkage = {0, 1, 1, 2, 2, 3};
  const std::vector<unsigned> side_set_flag = {0, 0};
  const std::vector<unsigned> side_node_count = {1, 1};
  const std::vector<unsigned> side_to_node_linkage = {0, 3};
  const std::vector<double> coordinates = {0.0, 1.0, 2.0, 3.0};
  const std::vector<unsigned> global_node_number = {0, 1, 2, 3};
  const std::vector<unsigned> num_nodes_per_face_per_cell = {1, 1, 1, 1, 1, 1};

  std::shared_ptr<rtt_mesh::Draco_Mesh> const global_mesh(new rtt_mesh::Draco_Mesh(
      dim, geom, num_faces_per_cell, cell_to_node_linkage, side_set_flag, side_node_count,
      side_to_node_linkage, coordinates, global_node_number, num_nodes_per_face_per_cell));

  // set global cell partition with 1st cell on rank 0 and the others on rank 1
  const std::vector<std::vector<int64_t>> global_cell_part = {{0}, {1, 2}};

  // create a mesh splitter object
  rtt_mesh::Draco_Mesh_Splitter mesh_splitter(global_cell_part);

  // create local mesh
  std::shared_ptr<rtt_mesh::Draco_Mesh> local_mesh = mesh_splitter.split_mesh(global_mesh);

  // check properties of child meshes are the same as parent
  FAIL_IF_NOT(local_mesh->get_dimension() == 1);
  FAIL_IF_NOT(local_mesh->get_geometry() == geom);
  FAIL_IF_NOT(local_mesh->get_ghost_cell_numbers().size() == 1);
  FAIL_IF_NOT(local_mesh->get_ghost_cell_numbers()[0] == 1);

  // check local mesh properties
  if (rtt_c4::node() == 0) {
    FAIL_IF_NOT(local_mesh->get_num_cells() == 1);
    FAIL_IF_NOT(local_mesh->get_num_nodes() == 2);
    FAIL_IF_NOT(local_mesh->get_ghost_cell_ranks()[0] == 1);
    FAIL_IF_NOT(local_mesh->get_global_nodes().size() == 2);
    FAIL_IF_NOT(local_mesh->get_global_nodes()[0] == 0);
    FAIL_IF_NOT(local_mesh->get_global_nodes()[1] == 1);
    FAIL_IF_NOT(local_mesh->get_node_coord_vec().size() == 2);
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(local_mesh->get_node_coord_vec()[0][0], 0.0));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(local_mesh->get_node_coord_vec()[1][0], 1.0));
  } else {
    FAIL_IF_NOT(local_mesh->get_num_cells() == 2);
    FAIL_IF_NOT(local_mesh->get_num_nodes() == 3);
    FAIL_IF_NOT(local_mesh->get_ghost_cell_ranks()[0] == 0);
    FAIL_IF_NOT(local_mesh->get_global_nodes().size() == 3);
    FAIL_IF_NOT(local_mesh->get_global_nodes()[0] == 1);
    FAIL_IF_NOT(local_mesh->get_global_nodes()[1] == 2);
    FAIL_IF_NOT(local_mesh->get_global_nodes()[2] == 3);
    FAIL_IF_NOT(local_mesh->get_node_coord_vec().size() == 3);
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(local_mesh->get_node_coord_vec()[0][0], 1.0));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(local_mesh->get_node_coord_vec()[1][0], 2.0));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(local_mesh->get_node_coord_vec()[2][0], 3.0));
  }

  // now apply a discontiguous global cell partition
  const std::vector<std::vector<int64_t>> global_cell_part2 = {{0, 2}, {1}};

  // reset the partition in the current mesh_splitter object
  mesh_splitter.reset_global_cells_per_rank(global_cell_part2);

  // recreate local mesh
  local_mesh = mesh_splitter.split_mesh(global_mesh);

  // check local mesh properties
  if (rtt_c4::node() == 0) {
    FAIL_IF_NOT(local_mesh->get_num_cells() == 2);
    FAIL_IF_NOT(local_mesh->get_num_nodes() == 4);
    FAIL_IF_NOT(local_mesh->get_ghost_cell_ranks().size() == 2);
    FAIL_IF_NOT(local_mesh->get_ghost_cell_ranks()[0] == 1);
    FAIL_IF_NOT(local_mesh->get_ghost_cell_ranks()[1] == 1);
    FAIL_IF_NOT(local_mesh->get_global_nodes().size() == 4);
    FAIL_IF_NOT(local_mesh->get_global_nodes()[0] == 0);
    FAIL_IF_NOT(local_mesh->get_global_nodes()[1] == 1);
    FAIL_IF_NOT(local_mesh->get_global_nodes()[2] == 2);
    FAIL_IF_NOT(local_mesh->get_global_nodes()[3] == 3);
    FAIL_IF_NOT(local_mesh->get_node_coord_vec().size() == 4);
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(local_mesh->get_node_coord_vec()[0][0], 0.0));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(local_mesh->get_node_coord_vec()[1][0], 1.0));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(local_mesh->get_node_coord_vec()[2][0], 2.0));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(local_mesh->get_node_coord_vec()[3][0], 3.0));
  } else {
    FAIL_IF_NOT(local_mesh->get_num_cells() == 1);
    FAIL_IF_NOT(local_mesh->get_num_nodes() == 2);
    FAIL_IF_NOT(local_mesh->get_ghost_cell_ranks().size() == 2);
    FAIL_IF_NOT(local_mesh->get_ghost_cell_ranks()[0] == 0);
    FAIL_IF_NOT(local_mesh->get_ghost_cell_ranks()[1] == 0);
    FAIL_IF_NOT(local_mesh->get_global_nodes().size() == 2);
    FAIL_IF_NOT(local_mesh->get_global_nodes()[0] == 1);
    FAIL_IF_NOT(local_mesh->get_global_nodes()[1] == 2);
    FAIL_IF_NOT(local_mesh->get_node_coord_vec().size() == 2);
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(local_mesh->get_node_coord_vec()[0][0], 1.0));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(local_mesh->get_node_coord_vec()[1][0], 2.0));
  }

  // now split a vector over global cell index into rank-local vectors
  const std::vector<double> global_scalars = {-2.0, 1.5,
                                              3.1416}; // NOLINT(modernize-use-std-numbers)

  // create local vector
  const auto local_scalars = mesh_splitter.split_vector(global_scalars);

  if (rtt_c4::node() == 0) {
    FAIL_IF_NOT(local_scalars.size() == 2);
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(local_scalars[0], -2.0));
    FAIL_IF_NOT(
        rtt_dsxx::soft_equiv(local_scalars[1], 3.1416)); // NOLINT(modernize-use-std-numbers)
  } else {
    FAIL_IF_NOT(local_scalars.size() == 1);
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(local_scalars[0], 1.5));
  }

  // successful test output
  if (ut.numFails == 0)
    PASSMSG("1D Draco_Mesh_Splitter tests ok.");
  return;
}

void split_3d_mesh_dd(rtt_c4::ParallelUnitTest &ut) {

  const std::string inputpath = ut.getTestSourcePath();
  const std::string filename = inputpath + "/x3d.mesh3d_8cell.in";
  const std::vector<std::string> boundary_filenames = {
      inputpath + "/x3d.mesh3d_8cell.bdy1.in", inputpath + "/x3d.mesh3d_8cell.bdy2.in",
      inputpath + "/x3d.mesh3d_8cell.bdy3.in", inputpath + "/x3d.mesh3d_8cell.bdy4.in",
      inputpath + "/x3d.mesh3d_8cell.bdy5.in", inputpath + "/x3d.mesh3d_8cell.bdy6.in"};

  // set some arbitrary boundary files
  const std::vector<unsigned> boundary_flags = {3, 1, 0, 2, 5, 4};

  // construct reader
  std::shared_ptr<rtt_mesh::X3D_Draco_Mesh_Reader> const x3d_reader(
      new rtt_mesh::X3D_Draco_Mesh_Reader(filename, boundary_filenames, boundary_flags));

  // read mesh
  x3d_reader->read_mesh();

  // planar geometry
  const rtt_mesh_element::Geometry geom = rtt_mesh_element::Geometry::CARTESIAN;

  // instantiate a mesh builder and build the mesh
  rtt_mesh::Draco_Mesh_Builder<rtt_mesh::X3D_Draco_Mesh_Reader> mesh_builder(x3d_reader);
  std::shared_ptr<rtt_mesh::Draco_Mesh> const global_mesh = mesh_builder.build_mesh(geom);

  // set a partition into 4 2x1x1 cell domains
  const std::vector<std::vector<int64_t>> global_cell_part = {{0, 1}, {2, 3}, {4, 5}, {6, 7}};

  // create a mesh splitter object
  rtt_mesh::Draco_Mesh_Splitter mesh_splitter(global_cell_part);

  // create local mesh
  std::shared_ptr<rtt_mesh::Draco_Mesh> local_mesh = mesh_splitter.split_mesh(global_mesh);

  // check properties common among ranks
  FAIL_IF_NOT(local_mesh->get_dimension() == 3);
  FAIL_IF_NOT(local_mesh->get_geometry() == geom);
  FAIL_IF_NOT(local_mesh->get_num_cells() == 2);
  FAIL_IF_NOT(local_mesh->get_num_nodes() == 12);
  FAIL_IF_NOT(local_mesh->get_ghost_cell_ranks().size() == 4);
  FAIL_IF_NOT(local_mesh->get_global_nodes().size() == 12);
  FAIL_IF_NOT(local_mesh->get_node_coord_vec().size() == 12);

  // check global node indices per rank
  if (rtt_c4::node() == 0) {
    const std::vector<unsigned> global_node_ref = {0, 1, 2, 3, 4, 5, 9, 10, 11, 12, 13, 14};
    for (int i = 0; i < 12; ++i)
      FAIL_IF_NOT(local_mesh->get_global_node_index(i) == global_node_ref[i]);
  } else if (rtt_c4::node() == 3) {
    const std::vector<unsigned> global_node_ref = {12, 13, 14, 15, 16, 17, 21, 22, 23, 24, 25, 26};
    for (int i = 0; i < 12; ++i)
      FAIL_IF_NOT(local_mesh->get_global_node_index(i) == global_node_ref[i]);
  }

  // now apply a discontiguous global cell partition
  const std::vector<std::vector<int64_t>> global_cell_part2 = {{0, 7}, {1, 6}, {2, 5}, {3, 4}};

  // reset the partition in the current mesh_splitter object
  mesh_splitter.reset_global_cells_per_rank(global_cell_part2);

  // recreate local mesh
  local_mesh = mesh_splitter.split_mesh(global_mesh);

  FAIL_IF_NOT(local_mesh->get_num_cells() == 2);
  FAIL_IF_NOT(local_mesh->get_num_nodes() == 15);
  FAIL_IF_NOT(local_mesh->get_ghost_cell_ranks().size() == 6);
  FAIL_IF_NOT(local_mesh->get_global_nodes().size() == 15);
  FAIL_IF_NOT(local_mesh->get_node_coord_vec().size() == 15);

  // check coordinates on a rank
  if (rtt_c4::node() == 1) {
    const std::vector<std::vector<double>> ref_coords = {
        {0.5, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.5, 1.0, 0.0}, {1.0, 1.0, 0.0}, {0.5, 0.0, 2.0},
        {1.0, 0.0, 2.0}, {0.0, 1.0, 2.0}, {0.5, 1.0, 2.0}, {1.0, 1.0, 2.0}, {0.0, 2.0, 2.0},
        {0.5, 2.0, 2.0}, {0.0, 1.0, 4.0}, {0.5, 1.0, 4.0}, {0.0, 2.0, 4.0}, {0.5, 2.0, 4.0}};
    const auto coords = local_mesh->get_node_coord_vec();
    for (int i = 0; i < 15; ++i) {
      for (int j = 0; j < 3; ++j) {
        FAIL_IF_NOT(rtt_dsxx::soft_equiv(coords[i][j], ref_coords[i][j]));
      }
    }
  }

  // now split a vector over global cell index into rank-local vectors
  const std::vector<double> global_scalars = {1.1, 2.2, 3.4, 4.5, 2.3, 3.5, 4.7, 5.001};

  // create local vector
  const auto local_scalars = mesh_splitter.split_vector(global_scalars);

  FAIL_IF_NOT(local_scalars.size() == 2);
  if (rtt_c4::node() == 0) {
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(local_scalars[0], 1.1));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(local_scalars[1], 5.001));
  } else if (rtt_c4::node() == 1) {
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(local_scalars[0], 2.2));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(local_scalars[1], 4.7));
  } else if (rtt_c4::node() == 2) {
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(local_scalars[0], 3.4));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(local_scalars[1], 3.5));
  } else if (rtt_c4::node() == 3) {
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(local_scalars[0], 4.5));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(local_scalars[1], 2.3));
  }

  // successful test output
  if (ut.numFails == 0)
    PASSMSG("3D Draco_Mesh_Splitter tests ok.");
  return;
}

} // namespace

//------------------------------------------------------------------------------------------------//

int main(int argc, char *argv[]) {
  rtt_c4::ParallelUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    if (rtt_c4::nodes() == 2) {
      split_1d_mesh_dd(ut);
      //split_2d_mesh_dd(ut);
    } else if (rtt_c4::nodes() == 4) {
      split_3d_mesh_dd(ut);
    } else {
      Insist(false, "This test only uses 2 or 4 PE.");
    }
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of mesh/test/tstDraco_Mesh_Splitter.cc
//------------------------------------------------------------------------------------------------//
