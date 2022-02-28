//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   mesh/test/tstX3D_Draco_Mesh_Reader.cc
 * \author Ryan Wollaeger <wollaeger@lanl.gov>
 * \date   Tuesday, Jul 10, 2018, 10:23 am
 * \brief  X3D_Draco_Mesh_Reader class unit test.
 * \note   Copyright (C) 2018-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "Test_Mesh_Interface.hh"
#include "c4/ParallelUnitTest.hh"
#include "ds++/Release.hh"
#include "ds++/Soft_Equivalence.hh"
#include "mesh/Draco_Mesh_Builder.hh"
#include "mesh/X3D_Draco_Mesh_Reader.hh"

using rtt_mesh::Draco_Mesh;
using rtt_mesh::Draco_Mesh_Builder;
using rtt_mesh::X3D_Draco_Mesh_Reader;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

// Parse an X3D file format and compare to reference
void read_x3d_mesh_2d(rtt_c4::ParallelUnitTest &ut) {

  // >>> PARSE MESH

  const std::string inputpath = ut.getTestSourcePath();
  const std::string filename = inputpath + "x3d.mesh.in";
  const std::vector<std::string> bdy_filenames = {
      inputpath + "x3d.mesh.bdy1.in", inputpath + "x3d.mesh.bdy2.in",
      inputpath + "x3d.mesh.bdy3.in", inputpath + "x3d.mesh.bdy4.in"};
  const std::vector<unsigned> bdy_flags = {3, 1, 0, 2};

  // construct reader
  std::shared_ptr<X3D_Draco_Mesh_Reader> x3d_reader(
      new X3D_Draco_Mesh_Reader(filename, bdy_filenames, bdy_flags));

  // read mesh
  x3d_reader->read_mesh();

  // >>> CHECK HEADER DATA

  FAIL_IF_NOT(x3d_reader->get_process() == 0);
  FAIL_IF_NOT(x3d_reader->get_numdim() == 2);
  FAIL_IF_NOT(x3d_reader->get_numcells() == 1);
  FAIL_IF_NOT(x3d_reader->get_numnodes() == 4);

  // >>> CHECK CELL-NODE DATA

  FAIL_IF_NOT(x3d_reader->get_celltype(0) == 4);

  std::vector<unsigned> test_cellnodes = {0, 1, 1, 3, 3, 2, 2, 0};
  FAIL_IF_NOT(x3d_reader->get_cellnodes(0) == test_cellnodes);

  // >>> CHECK NODE-COORD DATA

  std::vector<std::vector<double>> test_coords = {
      {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {1.0, 1.0, 0.0}};
  for (int node = 0; node < 4; ++node)
    FAIL_IF_NOT(x3d_reader->get_nodecoord(node) == test_coords[node]);

  // >>> CHECK SIDE DATA

  FAIL_IF_NOT(x3d_reader->get_numsides() == 4);

  std::vector<std::vector<unsigned>> test_sidenodes = {{0, 1}, {1, 3}, {3, 2}, {2, 0}};

  // check each side's data
  for (int side = 0; side < 4; ++side) {

    // sides must always give 2 nodes per face in X3D
    FAIL_IF_NOT(x3d_reader->get_sidetype(side) == 2);

    // boundary conditions are not supplied in X3D (note this check is specialized for the 1-cell
    // mesh)
    FAIL_IF_NOT(x3d_reader->get_sideflag(side) == bdy_flags[side]);

    // check node indices
    FAIL_IF_NOT(x3d_reader->get_sidenodes(side) == test_sidenodes[side]);
  }

  // >>> CHECK BC-NODE MAP

  const std::map<size_t, std::vector<unsigned>> &bc_node_map = x3d_reader->get_bc_node_map();

  FAIL_IF_NOT(bc_node_map.size() == 4);

  std::vector<std::vector<unsigned>> test_bc_nodes = {{0, 1}, {1, 3}, {2, 3}, {0, 2}};

  for (size_t ibc = 0; ibc < 4; ++ibc) {
    FAIL_IF_NOT(bc_node_map.at(ibc) == test_bc_nodes[ibc]);
  }

  // successful test output
  if (ut.numFails == 0)
    PASSMSG("2D X3D_Draco_Mesh_Reader parsing tests ok.");
  return;
}

//------------------------------------------------------------------------------------------------//
// Parse and build an X3D file format and compare to reference mesh
void build_x3d_mesh_2d(rtt_c4::ParallelUnitTest &ut) {

  // use Cartesian geometry
  const Draco_Mesh::Geometry geometry = Draco_Mesh::Geometry::CARTESIAN;

  // >>> CREATE REFERENCE MESH

  // set the number of cells and nodes
  const size_t num_xdir = 1;
  const size_t num_ydir = 1;

  // generate a constainer for data needed in mesh construction
  rtt_mesh_test::Test_Mesh_Interface mesh_iface(num_xdir, num_ydir, {}, 0.0, 0.0);

  // short-cut to some arrays
  const std::vector<unsigned> &cell_type = mesh_iface.cell_type;
  const std::vector<unsigned> &cell_to_node_linkage = mesh_iface.cell_to_node_linkage;
  const std::vector<unsigned> &side_node_count = mesh_iface.side_node_count;
  const std::vector<unsigned> &side_to_node_linkage = mesh_iface.side_to_node_linkage;

  // instantiate the mesh
  std::shared_ptr<Draco_Mesh> ref_mesh(
      new Draco_Mesh(mesh_iface.dim, geometry, cell_type, cell_to_node_linkage,
                     mesh_iface.side_set_flag, side_node_count, side_to_node_linkage,
                     mesh_iface.coordinates, mesh_iface.global_node_number, mesh_iface.face_type));

  // >>> PARSE AND BUILD MESH

  const std::string inputpath = ut.getTestSourcePath();
  const std::string filename = inputpath + "x3d.mesh.in";
  const std::vector<std::string> bdy_filenames = {
      inputpath + "x3d.mesh.bdy1.in", inputpath + "x3d.mesh.bdy2.in",
      inputpath + "x3d.mesh.bdy3.in", inputpath + "x3d.mesh.bdy4.in"};

  // construct reader
  std::shared_ptr<X3D_Draco_Mesh_Reader> x3d_reader(
      new X3D_Draco_Mesh_Reader(filename, bdy_filenames, {}));

  // read mesh
  x3d_reader->read_mesh();

  // instantiate a mesh builder and build the mesh
  Draco_Mesh_Builder<X3D_Draco_Mesh_Reader> mesh_builder(x3d_reader);
  std::shared_ptr<Draco_Mesh> mesh = mesh_builder.build_mesh(geometry);

  // check that the scalar data is correct
  FAIL_IF(mesh->get_dimension() != ref_mesh->get_dimension());
  FAIL_IF(mesh->get_geometry() != ref_mesh->get_geometry());
  FAIL_IF(mesh->get_num_cells() != ref_mesh->get_num_cells());
  FAIL_IF(mesh->get_num_nodes() != ref_mesh->get_num_nodes());

  // check that layout is correct (empty for one cell, no side or ghost data)
  FAIL_IF((mesh->get_cc_linkage()).size() > 0);
  FAIL_IF((mesh->get_cs_linkage()).size() != 1);
  FAIL_IF((mesh->get_cg_linkage()).size() > 0);

  // check side flag indices (should be different)
  FAIL_IF(mesh->get_side_set_flag() == ref_mesh->get_side_set_flag());

  // check ghost cell data (should be empty defaults)
  FAIL_IF(mesh->get_ghost_cell_numbers() != ref_mesh->get_ghost_cell_numbers());
  FAIL_IF(mesh->get_ghost_cell_ranks() != ref_mesh->get_ghost_cell_ranks());

  // check that the vector of coordinates match the reference mesh
  FAIL_IF(!rtt_dsxx::soft_equiv(mesh->get_node_coord_vec(), ref_mesh->get_node_coord_vec()));

  // check that each cell has the correct sides
  {
    std::vector<unsigned> test_sn_linkage = mesh_iface.flatten_sn_linkage(mesh->get_cs_linkage());

    // check that sn_linkage is a permutation of the original side-node linkage
    auto sn_first = side_to_node_linkage.begin();
    auto test_sn_first = test_sn_linkage.begin();

    for (unsigned side = 0; side < mesh_iface.num_sides; ++side) {

      // check that sn_linkage is a permutation of the original side-node linkage
      FAIL_IF(!std::is_permutation(test_sn_first, test_sn_first + side_node_count[side], sn_first,
                                   sn_first + side_node_count[side]));

      // update the iterators
      sn_first += side_node_count[side];
      test_sn_first += side_node_count[side];
    }
  }

  // successful test output
  if (ut.numFails == 0)
    PASSMSG("2D X3D_Draco_Mesh_Reader mesh build tests ok.");
  return;
}

//----------------------------------------------------------------------------------------//
// Parse a larger Voronoi X3D file format and compare to reference
void read_voronoi_mesh(rtt_c4::ParallelUnitTest &ut) {

  // >>> PARSE MESH

  const std::string inputpath = ut.getTestSourcePath();
  const std::string filename = inputpath + "pipe_voronoi.x3d.in";
  const std::vector<std::string> bdy_filenames;
  const std::vector<unsigned> bdy_flags;

  // construct reader
  std::shared_ptr<X3D_Draco_Mesh_Reader> x3d_reader(
      new X3D_Draco_Mesh_Reader(filename, bdy_filenames, bdy_flags));

  // read mesh
  x3d_reader->read_mesh();

  // >>> CHECK HEADER DATA

  FAIL_IF_NOT(x3d_reader->get_process() == 0);
  FAIL_IF_NOT(x3d_reader->get_numdim() == 2);
  FAIL_IF_NOT(x3d_reader->get_numcells() == 1590);
  FAIL_IF_NOT(x3d_reader->get_numnodes() == 2291);

  // >>> CHECK CELL-NODE DATA

  FAIL_IF_NOT(x3d_reader->get_celltype(0) == 5);

  std::vector<unsigned> test_cellnodes = {0, 938, 938, 212, 212, 199, 199, 939, 939, 0};
  FAIL_IF_NOT(x3d_reader->get_cellnodes(0) == test_cellnodes);

  // CHECK MATIDS
  FAIL_IF_NOT(x3d_reader->get_matid(0) == "1");

  // successful test output
  if (ut.numFails == 0)
    PASSMSG("2D X3D_Draco_Mesh_Reader parsing tests ok.");
  return;
}

//------------------------------------------------------------------------------------------------//
// Parse an X3D file format and compare to reference
void read_x3d_mesh_3d(rtt_c4::ParallelUnitTest &ut) {

  // >>> PARSE MESH

  const std::string inputpath = ut.getTestSourcePath();
  const std::string filename = inputpath + "x3d.mesh3d_1cell.in";
  const std::vector<std::string> bdy_filenames = {
      inputpath + "x3d.mesh3d_1cell.bdy1.in", inputpath + "x3d.mesh3d_1cell.bdy2.in",
      inputpath + "x3d.mesh3d_1cell.bdy3.in", inputpath + "x3d.mesh3d_1cell.bdy4.in",
      inputpath + "x3d.mesh3d_1cell.bdy5.in", inputpath + "x3d.mesh3d_1cell.bdy6.in"};
  const std::vector<unsigned> bdy_flags = {3, 1, 0, 2, 1, 5};

  // construct reader
  std::shared_ptr<X3D_Draco_Mesh_Reader> x3d_reader(
      new X3D_Draco_Mesh_Reader(filename, bdy_filenames, bdy_flags));

  // read mesh
  x3d_reader->read_mesh();

  // >>> CHECK HEADER DATA

  FAIL_IF_NOT(x3d_reader->get_process() == 0);
  FAIL_IF_NOT(x3d_reader->get_numdim() == 3);
  FAIL_IF_NOT(x3d_reader->get_numcells() == 1);
  FAIL_IF_NOT(x3d_reader->get_numnodes() == 8);

  // >>> CHECK CELL-NODE DATA

  FAIL_IF_NOT(x3d_reader->get_celltype(0) == 6);

  std::vector<unsigned> test_cellnodes = {0, 2, 3, 1, 5, 7, 6, 4, 0, 1, 5, 4,
                                          6, 7, 3, 2, 0, 4, 6, 2, 3, 7, 5, 1};
  FAIL_IF_NOT(x3d_reader->get_cellnodes(0) == test_cellnodes);

  // >>> CHECK NODE-COORD DATA

  std::vector<std::vector<double>> test_coords = {{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0},
                                                  {1.0, 1.0, 0.0}, {0.0, 0.0, 1.0}, {1.0, 0.0, 1.0},
                                                  {0.0, 1.0, 1.0}, {1.0, 1.0, 1.0}};
  for (int node = 0; node < 8; ++node)
    FAIL_IF_NOT(x3d_reader->get_nodecoord(node) == test_coords[node]);

  // >>> CHECK SIDE DATA

  FAIL_IF_NOT(x3d_reader->get_numsides() == 6);

  // set test boundary nodes
  std::vector<std::vector<unsigned>> test_sidenodes = {{0, 4, 6, 2}, {3, 7, 5, 1}, {0, 1, 5, 4},
                                                       {6, 7, 3, 2}, {0, 2, 3, 1}, {5, 7, 6, 4}};

  // check each side's data
  for (int side = 0; side < 6; ++side) {

    // sides must always give 2 nodes per face in X3D
    FAIL_IF_NOT(x3d_reader->get_sidetype(side) == 4);

    // boundary conditions are not supplied in X3D (note this check is specialized for the 1-cell
    // mesh)
    FAIL_IF_NOT(x3d_reader->get_sideflag(side) == bdy_flags[side]);

    // check node indices
    FAIL_IF_NOT(x3d_reader->get_sidenodes(side) == test_sidenodes[side]);
  }

  // >>> CHECK BC-NODE MAP

  const std::map<size_t, std::vector<unsigned>> &bc_node_map = x3d_reader->get_bc_node_map();

  FAIL_IF_NOT(bc_node_map.size() == 6);

  std::vector<std::vector<unsigned>> test_bc_nodes = {{0, 2, 4, 6}, {1, 3, 5, 7}, {0, 1, 4, 5},
                                                      {2, 3, 6, 7}, {0, 1, 2, 3}, {4, 5, 6, 7}};

  for (size_t ibc = 0; ibc < 6; ++ibc) {
    FAIL_IF_NOT(bc_node_map.at(ibc) == test_bc_nodes[ibc]);
  }

  // successful test output
  if (ut.numFails == 0)
    PASSMSG("3D X3D_Draco_Mesh_Reader parsing tests ok.");
  return;
}

//------------------------------------------------------------------------------------------------//
// Parse and build an X3D file format and compare to reference mesh
void build_x3d_mesh_3d(rtt_c4::ParallelUnitTest &ut) {

  // use Cartesian geometry
  const Draco_Mesh::Geometry geometry = Draco_Mesh::Geometry::CARTESIAN;

  // >>> PARSE MESH

  const std::string inputpath = ut.getTestSourcePath();
  const std::string filename = inputpath + "x3d.mesh3d_1cell.in";
  const std::vector<std::string> bdy_filenames = {
      inputpath + "x3d.mesh3d_1cell.bdy1.in", inputpath + "x3d.mesh3d_1cell.bdy2.in",
      inputpath + "x3d.mesh3d_1cell.bdy3.in", inputpath + "x3d.mesh3d_1cell.bdy4.in",
      inputpath + "x3d.mesh3d_1cell.bdy5.in", inputpath + "x3d.mesh3d_1cell.bdy6.in"};
  const std::vector<unsigned> bdy_flags = {3, 1, 0, 2, 1, 5};

  // construct reader
  std::shared_ptr<X3D_Draco_Mesh_Reader> x3d_reader(
      new X3D_Draco_Mesh_Reader(filename, bdy_filenames, bdy_flags));

  // read mesh
  x3d_reader->read_mesh();

  // instantiate a mesh builder and build the mesh
  Draco_Mesh_Builder<X3D_Draco_Mesh_Reader> mesh_builder(x3d_reader);
  std::shared_ptr<Draco_Mesh> mesh = mesh_builder.build_mesh(geometry);

  // check that the scalar data is correct
  FAIL_IF(mesh->get_dimension() != 3);
  FAIL_IF(mesh->get_geometry() != Draco_Mesh::Geometry::CARTESIAN);
  FAIL_IF(mesh->get_num_cells() != 1);
  FAIL_IF(mesh->get_num_nodes() != 8);

  // check layout sizes
  FAIL_IF((mesh->get_cc_linkage()).size() > 0);
  FAIL_IF((mesh->get_cs_linkage()).size() != 1);
  FAIL_IF((mesh->get_cg_linkage()).size() > 0);

  // get the boundary cell layout
  Draco_Mesh::Layout cellside_linkage = mesh->get_cs_linkage();

  // check that 6 faces are adjacent to the side
  FAIL_IF_NOT(cellside_linkage.at(0).size() == 6);

  // test nodes ordered clockwise out of cell
  std::vector<std::vector<unsigned>> ref_nodes = {{0, 2, 3, 1}, {5, 7, 6, 4}, {0, 1, 5, 4},
                                                  {6, 7, 3, 2}, {0, 4, 6, 2}, {3, 7, 5, 1}};

  // test "side" index of face (side index ordered low-x, hig-x, low-y, hig-y, low-z, hig-z)
  // side indices are in order of parsered file, not necessarily same as face ordering in layout
  std::vector<unsigned> ref_sides = {4, 5, 2, 3, 0, 1};

  if (cellside_linkage.at(0).size() == 6) {

    // check number of nodes per face
    for (int face = 0; face < 6; ++face) {
      FAIL_IF_NOT(cellside_linkage.at(0)[face].first == ref_sides[face]);
      FAIL_IF_NOT(cellside_linkage.at(0)[face].second.size() == 4);
      FAIL_IF_NOT(cellside_linkage.at(0)[face].second == ref_nodes[face]);
    }
  }

  // successful test output
  if (ut.numFails == 0)
    PASSMSG("3D X3D_Draco_Mesh_Reader mesh build tests ok.");
  return;
}

//------------------------------------------------------------------------------------------------//
// Parse and build a 3D random mesh
void build_x3d_rnd_mesh_3d(rtt_c4::ParallelUnitTest &ut) {

  // use Cartesian geometry
  const Draco_Mesh::Geometry geometry = Draco_Mesh::Geometry::CARTESIAN;

  // >>> PARSE MESH

  const std::string inputpath = ut.getTestSourcePath();
  const std::string filename = inputpath + "x3d.rnd_mesh3d.in";
  const std::vector<std::string> bdy_filenames = {
      inputpath + "x3d.rnd_mesh3d.bdy1.in", inputpath + "x3d.rnd_mesh3d.bdy2.in",
      inputpath + "x3d.rnd_mesh3d.bdy3.in", inputpath + "x3d.rnd_mesh3d.bdy4.in",
      inputpath + "x3d.rnd_mesh3d.bdy5.in", inputpath + "x3d.rnd_mesh3d.bdy6.in"};
  const std::vector<unsigned> bdy_flags = {3, 1, 0, 2, 1, 5};

  // construct reader
  std::shared_ptr<X3D_Draco_Mesh_Reader> x3d_reader(
      new X3D_Draco_Mesh_Reader(filename, bdy_filenames, bdy_flags));

  // read mesh
  x3d_reader->read_mesh();

  // instantiate a mesh builder and build the mesh
  Draco_Mesh_Builder<X3D_Draco_Mesh_Reader> mesh_builder(x3d_reader);
  std::shared_ptr<Draco_Mesh> mesh = mesh_builder.build_mesh(geometry);

  // check that the scalar data is correct
  FAIL_IF(mesh->get_dimension() != 3);
  FAIL_IF(mesh->get_geometry() != Draco_Mesh::Geometry::CARTESIAN);
  FAIL_IF(mesh->get_num_cells() != 64);
  FAIL_IF(mesh->get_num_nodes() != 125);

  // check layout sizes
  FAIL_IF((mesh->get_cc_linkage()).size() != 64);
  FAIL_IF((mesh->get_cs_linkage()).size() != 56);
  FAIL_IF((mesh->get_cg_linkage()).size() > 0);

  // \todo: add more testing here

  // successful test output
  if (ut.numFails == 0)
    PASSMSG("3D X3D_Draco_Mesh_Reader random mesh build tests ok.");
  return;
}

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  rtt_c4::ParallelUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    Insist(rtt_c4::nodes() == 1, "This test only uses 1 PE.");
    read_x3d_mesh_2d(ut);
    build_x3d_mesh_2d(ut);
    read_voronoi_mesh(ut);
    read_x3d_mesh_3d(ut);
    build_x3d_mesh_3d(ut);
    build_x3d_rnd_mesh_3d(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of mesh/test/tstX3D_Draco_Mesh_Reader.cc
//------------------------------------------------------------------------------------------------//
