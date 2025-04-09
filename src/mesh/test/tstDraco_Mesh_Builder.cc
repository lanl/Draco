//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   mesh/test/tstDraco_Mesh_Builder.cc
 * \author Ryan Wollaeger <wollaeger@lanl.gov>
 * \date   Sunday, Jul 01, 2018, 18:21 pm
 * \brief  Draco_Mesh_Builder class unit test.
 * \note   Copyright (C) 2018-2025 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "Test_Mesh_Interface.hh"
#include "dsxx/Assert.hh"
#include "dsxx/Release.hh"
#include "dsxx/Soft_Equivalence.hh"
#include "dsxx/UnitTest.hh"
#include "c4/C4_Functions.hh"
#include "c4/ParallelUnitTest.hh"
#include "mesh/Draco_Mesh_Builder.hh"
#include "mesh/RTT_Draco_Mesh_Reader.hh"
#include "mesh_element/Geometry.hh"
#include <algorithm>
#include <cstddef>
#include <memory>
#include <string>
#include <vector>

using rtt_mesh::Draco_Mesh;
using rtt_mesh::Draco_Mesh_Builder;
using rtt_mesh::RTT_Draco_Mesh_Reader;

// cSpell:ignore xdir, ydir

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

// Parse a 2D Cartesian mesh and compare to reference
void build_cartesian_mesh_2d(rtt_c4::ParallelUnitTest &ut) {

  // use Cartesian geometry
  const Draco_Mesh::Geometry geometry = Draco_Mesh::Geometry::CARTESIAN;

  // >>> CREATE REFERENCE MESH

  // set the number of cells and nodes
  const size_t num_xdir = 2;
  const size_t num_ydir = 1;

  // generate a container for data needed in mesh construction
  rtt_mesh_test::Test_Mesh_Interface const mesh_iface(num_xdir, num_ydir);

  // short-cut to some arrays
  const std::vector<unsigned> &cell_type = mesh_iface.cell_type;
  const std::vector<unsigned> &face_type = mesh_iface.face_type;
  const std::vector<unsigned> &cell_to_node_linkage = mesh_iface.cell_to_node_linkage;
  const std::vector<unsigned> &side_node_count = mesh_iface.side_node_count;
  const std::vector<unsigned> &side_to_node_linkage = mesh_iface.side_to_node_linkage;

  // instantiate the mesh
  std::shared_ptr<Draco_Mesh> const ref_mesh(
      new Draco_Mesh(mesh_iface.dim, geometry, cell_type, cell_to_node_linkage,
                     mesh_iface.side_set_flag, side_node_count, side_to_node_linkage,
                     mesh_iface.coordinates, mesh_iface.global_node_number, face_type));

  // >>> PARSE AND BUILD MESH

  std::string const inputpath = ut.getTestSourcePath();
  std::string const filename = inputpath + "rttquad_2cell.mesh.in";

  // create an rtt mesh reader and read the mesh
  std::shared_ptr<RTT_Draco_Mesh_Reader> const rtt_mesh(new RTT_Draco_Mesh_Reader(filename));
  rtt_mesh->read_mesh();

  // instantiate a mesh builder and build the mesh
  Draco_Mesh_Builder<RTT_Draco_Mesh_Reader> mesh_builder(rtt_mesh);
  std::shared_ptr<Draco_Mesh> const mesh = mesh_builder.build_mesh(geometry);

  // check that the scalar data is correct
  FAIL_IF(mesh->get_dimension() != ref_mesh->get_dimension());
  FAIL_IF(mesh->get_geometry() != ref_mesh->get_geometry());
  FAIL_IF(mesh->get_num_cells() != ref_mesh->get_num_cells());
  FAIL_IF(mesh->get_num_nodes() != ref_mesh->get_num_nodes());

  // check side flag indices
  FAIL_IF(mesh->get_side_set_flag() != ref_mesh->get_side_set_flag());

  // check ghost cell data (should be empty defaults)
  FAIL_IF(mesh->get_ghost_cell_numbers() != ref_mesh->get_ghost_cell_numbers());
  FAIL_IF(mesh->get_ghost_cell_ranks() != ref_mesh->get_ghost_cell_ranks());

  // check that the vector of coordinates match the reference mesh
  FAIL_IF(!rtt_dsxx::soft_equiv(mesh->get_node_coord_vec(), ref_mesh->get_node_coord_vec()));

  // check that cell-to-node linkage data is correct
  {
    std::vector<unsigned> ref_cn_linkage = mesh_iface.flatten_cn_linkage(
        ref_mesh->get_cc_linkage(), ref_mesh->get_cs_linkage(), ref_mesh->get_cg_linkage());
    std::vector<unsigned> cn_linkage = mesh_iface.flatten_cn_linkage(
        mesh->get_cc_linkage(), mesh->get_cs_linkage(), mesh->get_cg_linkage());

    // check that cn_linkage is a permutation of the ref cell-node linkage
    auto ref_cn_first = ref_cn_linkage.begin();
    auto cn_first = cn_linkage.begin();

    unsigned fc_indx = 0;
    for (unsigned cell = 0; cell < mesh_iface.num_cells; ++cell) {

      unsigned num_nodes_per_cell = 0;
      for (unsigned face = 0; face < cell_type[cell]; ++face) {

        // accumulate number of nodes
        num_nodes_per_cell += face_type[fc_indx];

        // increment face-type counter
        fc_indx++;
      }

      // nodes must only be permuted at the cell level
      FAIL_IF_NOT(std::is_permutation(cn_first, cn_first + num_nodes_per_cell, ref_cn_first,
                                      ref_cn_first + num_nodes_per_cell));

      // update the iterators
      ref_cn_first += num_nodes_per_cell;
      cn_first += num_nodes_per_cell;
    }
  }

  // check that each cell has the correct sides
  {
    std::vector<unsigned> ref_sn_linkage =
        mesh_iface.flatten_sn_linkage(ref_mesh->get_cs_linkage());
    std::vector<unsigned> sn_linkage = mesh_iface.flatten_sn_linkage(mesh->get_cs_linkage());

    // check that sn_linkage is a permutation of the original side-node linkage
    auto ref_sn_first = ref_sn_linkage.begin();
    auto sn_first = sn_linkage.begin();
    for (unsigned side = 0; side < mesh_iface.num_sides; ++side) {

      // check that sn_linkage is a permutation of the original side-node linkage
      FAIL_IF_NOT(std::is_permutation(sn_first, sn_first + side_node_count[side], ref_sn_first,
                                      ref_sn_first + side_node_count[side]));

      // update the iterators
      ref_sn_first += side_node_count[side];
      sn_first += side_node_count[side];
    }
  }

  // build a modified mesh with new geometry, B.C.s, and coordinates
  {
    const Draco_Mesh::Geometry new_geometry = Draco_Mesh::Geometry::AXISYMMETRIC;

    // scale coordinates by constant scalar factor and reset them in mesh_builder
    constexpr double s = 3.0;
    const std::vector<double> coords = mesh_builder.get_coordinates();
    std::vector<double> scoords(coords.size());
    std::transform(coords.begin(), coords.end(), scoords.begin(), [](double v) { return s * v; });
    mesh_builder.reset_coordinates(scoords);
    FAIL_IF(rtt_dsxx::soft_equiv(coords, mesh_builder.get_coordinates()));

    // add 1 to B.C. side flags and reset them in mesh_builder
    constexpr unsigned j = 1;
    const std::vector<unsigned> bcs = mesh_builder.get_side_set_flag();
    std::vector<unsigned> bcsp1(bcs.size());
    std::transform(bcs.begin(), bcs.end(), bcsp1.begin(), [](unsigned v) { return v + j; });
    mesh_builder.reset_side_set_flag(bcsp1);
    FAIL_IF(bcs == mesh_builder.get_side_set_flag());

    // create a new mesh based on the old one, but with the above modifications
    std::shared_ptr<Draco_Mesh> const remesh = mesh_builder.rebuild_mesh(new_geometry);
    FAIL_IF(remesh->get_dimension() != ref_mesh->get_dimension());
    FAIL_IF(remesh->get_num_cells() != ref_mesh->get_num_cells());
    FAIL_IF(remesh->get_num_nodes() != ref_mesh->get_num_nodes());

    // check the geometry, coordinates, and B.C.s are modified
    FAIL_IF(remesh->get_geometry() == ref_mesh->get_geometry());
    FAIL_IF(remesh->get_side_set_flag() == ref_mesh->get_side_set_flag());
    FAIL_IF(remesh->get_node_coord_vec() == ref_mesh->get_node_coord_vec());
  }

  // successful test output
  if (ut.numFails == 0)
    PASSMSG("2D Cartesian Draco_Mesh_Builder tests okay.");
  return;
}

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  rtt_c4::ParallelUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    Insist(rtt_c4::nodes() == 1, "This test only uses 1 PE.");
    build_cartesian_mesh_2d(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of mesh/test/tstDraco_Mesh_Builder.cc
//------------------------------------------------------------------------------------------------//
