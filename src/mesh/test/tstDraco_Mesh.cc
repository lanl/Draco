//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   mesh/test/tstDraco_Mesh.cc
 * \author Ryan Wollaeger <wollaeger@lanl.gov>
 * \date   Thursday, Jun 07, 2018, 15:43 pm
 * \brief  Draco_Mesh class unit test.
 * \note   Copyright (C) 2018-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "Test_Mesh_Interface.hh"
#include "c4/ParallelUnitTest.hh"
#include "ds++/Release.hh"

using rtt_mesh::Draco_Mesh;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

// 1D Spherical mesh construction test
void spherical_mesh_1d(rtt_c4::ParallelUnitTest &ut) {

  // use Spherical geometry and 1D
  const unsigned dimension = 1;
  const Draco_Mesh::Geometry geometry = Draco_Mesh::Geometry::SPHERICAL;

  //>>> SET UP CELL AND NODE DATA

  // mesh:  |   |   |   |
  //       0.0 1.0 2.0 3.0

  const size_t num_cells = 3;
  const std::vector<unsigned> cell_type = {2, 2, 2};
  const std::vector<unsigned> cell_to_node_linkage = {0, 1, 1, 2, 2, 3};
  const std::vector<unsigned> side_set_flag = {0, 0};
  const std::vector<unsigned> side_node_count = {1, 1};
  const std::vector<unsigned> side_to_node_linkage = {0, 3};
  const std::vector<double> coordinates = {0.0, 1.0, 2.0, 3.0};
  const std::vector<unsigned> global_node_number = {0, 1, 2, 3};
  const std::vector<unsigned> face_type = {1, 1, 1, 1, 1, 1};

  // instantiate a mesh with face type input
  std::shared_ptr<Draco_Mesh> mesh(new Draco_Mesh(
      dimension, geometry, cell_type, cell_to_node_linkage, side_set_flag, side_node_count,
      side_to_node_linkage, coordinates, global_node_number, face_type));

  // instantiate a mesh without boundary faces supplied
  std::shared_ptr<Draco_Mesh> mesh_nobdy(
      new Draco_Mesh(dimension, geometry, cell_type, cell_to_node_linkage, {}, {}, {}, coordinates,
                     global_node_number, face_type));

  //>>> TEST THE MESH LAYOUTS

  // check that no ghost data has been generated
  FAIL_IF_NOT(mesh->get_cg_linkage().size() == 0);
  FAIL_IF_NOT(mesh_nobdy->get_cg_linkage().size() == 0);

  // get the cell-side layout generated by the mesh
  const Draco_Mesh::Layout bd_layout = mesh->get_cs_linkage();
  FAIL_IF_NOT(bd_layout.size() == 2);
  FAIL_IF_NOT(bd_layout.at(0).size() == 1);
  FAIL_IF_NOT(bd_layout.at(0)[0].first == 0);
  FAIL_IF_NOT(bd_layout.at(num_cells - 1).size() == 1);
  FAIL_IF_NOT(bd_layout.at(num_cells - 1)[0].first == 1);

  // check that the default boundary data was correctly generated
  const Draco_Mesh::Layout bd_layout_nobdy = mesh_nobdy->get_cs_linkage();
  FAIL_IF_NOT(bd_layout_nobdy == bd_layout);

  // get the cell-cell layout generated by the mesh
  const Draco_Mesh::Layout layout = mesh->get_cc_linkage();
  FAIL_IF_NOT(layout.size() == num_cells);

  // check inner boundary cell
  FAIL_IF_NOT(layout.at(0).size() == 1);
  FAIL_IF_NOT(layout.at(0)[0].first == 1);
  FAIL_IF_NOT(layout.at(0)[0].second.size() == 1);
  FAIL_IF_NOT(layout.at(0)[0].second[0] == 1);

  // check outer boundary cell
  FAIL_IF_NOT(layout.at(num_cells - 1).size() == 1);
  FAIL_IF_NOT(layout.at(num_cells - 1)[0].first == 1);
  FAIL_IF_NOT(layout.at(num_cells - 1)[0].second.size() == 1);
  FAIL_IF_NOT(layout.at(num_cells - 1)[0].second[0] == 2);

  // check internal connectivity
  for (unsigned cell = 1; cell < num_cells - 1; ++cell) {

    // check neighbor cells
    FAIL_IF_NOT(layout.at(cell).size() == 2);
    FAIL_IF_NOT(layout.at(cell)[0].first == cell - 1);
    FAIL_IF_NOT(layout.at(cell)[1].first == cell + 1);

    // check connecting nodes
    FAIL_IF_NOT(layout.at(cell)[0].second.size() == 1);
    FAIL_IF_NOT(layout.at(cell)[0].second[0] == cell_to_node_linkage[2 * cell]);
    FAIL_IF_NOT(layout.at(cell)[1].second.size() == 1);
    FAIL_IF_NOT(layout.at(cell)[1].second[0] == cell_to_node_linkage[2 * cell + 1]);

    // check face index of neighboring cell
    if (cell == 1) {
      FAIL_IF_NOT(mesh->next_face(cell + 1, 1) == 1);
    } else {
      FAIL_IF_NOT(mesh->next_face(cell + 1, 1) == 2);
      FAIL_IF_NOT(mesh->next_face(cell + 1, 2) == 1);
    }
  }

  // get the layout generated by the mesh
  const Draco_Mesh::Layout layout_nobdy = mesh->get_cc_linkage();
  FAIL_IF_NOT(layout_nobdy == layout);

  // successful test output
  if (ut.numFails == 0)
    PASSMSG("1D Draco_Mesh tests okay.");
  return;
}

// 2D Cartesian mesh construction test
void cartesian_mesh_2d(rtt_c4::ParallelUnitTest &ut) {

  // use Cartesian geometry
  const Draco_Mesh::Geometry geometry = Draco_Mesh::Geometry::CARTESIAN;

  //>>> SET UP CELL AND NODE DATA

  // set the number of cells and nodes
  const size_t num_xdir = 6;
  const size_t num_ydir = 5;

  // generate a container for data needed in mesh construction
  rtt_mesh_test::Test_Mesh_Interface mesh_iface(num_xdir, num_ydir);

  // short-cut to some arrays
  const std::vector<unsigned> &cell_type = mesh_iface.cell_type;
  const std::vector<unsigned> &cell_to_node_linkage = mesh_iface.cell_to_node_linkage;
  const std::vector<unsigned> &side_node_count = mesh_iface.side_node_count;
  const std::vector<unsigned> &side_to_node_linkage = mesh_iface.side_to_node_linkage;

  // instantiate a mesh with face type input
  std::shared_ptr<Draco_Mesh> mesh(new Draco_Mesh(
      mesh_iface.dim, geometry, mesh_iface.cell_type, mesh_iface.cell_to_node_linkage,
      mesh_iface.side_set_flag, mesh_iface.side_node_count, mesh_iface.side_to_node_linkage,
      mesh_iface.coordinates, mesh_iface.global_node_number, mesh_iface.face_type));

  // check that the scalar data is correct
  FAIL_IF_NOT(mesh->get_dimension() == 2);
  FAIL_IF_NOT(mesh->get_geometry() == Draco_Mesh::Geometry::CARTESIAN);
  FAIL_IF_NOT(mesh->get_num_cells() == mesh_iface.num_cells);
  FAIL_IF_NOT(mesh->get_num_nodes() == mesh_iface.num_nodes);

  // check that flat cell type and cell to node linkage is correct
  FAIL_IF_NOT(mesh->get_num_faces_per_cell() == cell_type);
  FAIL_IF_NOT(mesh->get_flat_cell_node_linkage() == cell_to_node_linkage);

  // check that flat side type and side to node linkage is correct
  FAIL_IF_NOT(mesh->get_side_node_count() == side_node_count);
  FAIL_IF_NOT(mesh->get_side_to_node_linkage() == side_to_node_linkage);
  FAIL_IF_NOT(mesh->get_side_set_flag() == mesh_iface.side_set_flag);

  // get the layout generated by the mesh
  const Draco_Mesh::Layout layout = mesh->get_cc_linkage();

  // check that the layout has been generated
  FAIL_IF_NOT(layout.size() == mesh_iface.num_cells);

  // check that each cell has the correct neighbors
  {
    std::map<unsigned, std::vector<unsigned>> test_cell_map;
    for (size_t j = 0; j < num_ydir; ++j) {
      for (size_t i = 0; i < num_xdir; ++i) {

        // calculate the cell index
        auto cell = static_cast<unsigned>(i + j * num_xdir);

        // calculate neighbor cell indices
        if (j > 0)
          test_cell_map[cell].push_back(cell - static_cast<unsigned>(num_xdir));
        if (i < num_xdir - 1)
          test_cell_map[cell].push_back(cell + 1);
        if (j < num_ydir - 1)
          test_cell_map[cell].push_back(cell + static_cast<unsigned>(num_xdir));
        if (i > 0)
          test_cell_map[cell].push_back(cell - 1);
      }
    }

    for (unsigned cell = 0; cell < mesh_iface.num_cells; ++cell) {

      // get number of faces per cell in layout
      const size_t num_faces = layout.at(cell).size();

      // check that the number of faces per cell is correct
      FAIL_IF_NOT(num_faces == test_cell_map[cell].size());

      // check that cell neighbors are correct
      for (unsigned face = 0; face < num_faces; ++face)
        FAIL_IF_NOT(layout.at(cell)[face].first == test_cell_map[cell][face]);
    }
  }

  // get the boundary layout generated by the mesh
  const Draco_Mesh::Layout bd_layout = mesh->get_cs_linkage();

  // check that the boundary (or side) layout has been generated
  FAIL_IF_NOT(bd_layout.size() <= mesh_iface.num_cells);

  // get the ghost-cell layout generated by the mesh
  const Draco_Mesh::Layout go_layout = mesh->get_cg_linkage();

  // check that there are no ghost cells for this mesh
  FAIL_IF_NOT(go_layout.size() == 0);

  // check that cell-to-node linkage data is correct
  {
    std::vector<unsigned> test_cn_linkage =
        mesh_iface.flatten_cn_linkage(layout, bd_layout, go_layout);

    // check that cn_linkage is a permutation of the original cell-node linkage
    auto cn_first = cell_to_node_linkage.begin();
    auto test_cn_first = test_cn_linkage.begin();
    for (unsigned cell = 0; cell < mesh_iface.num_cells; ++cell) {

      // get the unique cell nodes
      std::vector<unsigned> cell_nodes = mesh->get_cell_nodes(cell);

      // assume the duplicate nodes multiply the stride by 2 (only true for 2D)
      const size_t nnpc = 2 * cell_nodes.size();

      // nodes must only be permuted at the cell level (assumes 2D for nnpc)
      FAIL_IF_NOT(
          std::is_permutation(test_cn_first, test_cn_first + nnpc, cn_first, cn_first + nnpc));

      // check that unique node entries from mesh are correct
      std::sort(cell_nodes.begin(), cell_nodes.end());
      std::vector<unsigned> cn_vec(cn_first, cn_first + nnpc);
      std::sort(cn_vec.begin(), cn_vec.end());
      auto last = std::unique(cn_vec.begin(), cn_vec.end());
      cn_vec.erase(last, cn_vec.end());
      FAIL_IF_NOT(cell_nodes == cn_vec);

      // update the iterators
      cn_first += nnpc;
      test_cn_first += nnpc;
    }

    FAIL_IF_NOT(cn_first == cell_to_node_linkage.end());
    FAIL_IF_NOT(test_cn_first == test_cn_linkage.end());
  }

  // check that each cell has the correct sides
  {
    std::vector<unsigned> test_sn_linkage = mesh_iface.flatten_sn_linkage(bd_layout);

    // check that sn_linkage is a permutation of the original side-node linkage
    auto sn_first = side_to_node_linkage.begin();
    auto test_sn_first = test_sn_linkage.begin();
    for (unsigned side = 0; side < mesh_iface.num_sides; ++side) {

      // sn_linkage must be a permutation of the original side-node linkage
      FAIL_IF_NOT(std::is_permutation(test_sn_first, test_sn_first + side_node_count[side],
                                      sn_first, sn_first + side_node_count[side]));

      // update the iterators
      sn_first += side_node_count[side];
      test_sn_first += side_node_count[side];
    }
  }

  // test default side (boundary) data
  {
    // set reference default side-set flags
    std::vector<unsigned> nobc_side_set_flag(mesh_iface.num_sides, 0);

    // instantiate a version of the mesh without side (b.c.) data
    std::shared_ptr<Draco_Mesh> mesh_no_bc_data(new Draco_Mesh(
        mesh_iface.dim, geometry, cell_type, cell_to_node_linkage, {}, {}, {},
        mesh_iface.coordinates, mesh_iface.global_node_number, mesh_iface.face_type));

    // check that default data has been initialized at mesh boundaries
    FAIL_IF_NOT(mesh_no_bc_data->get_side_node_count() == side_node_count);
    const std::vector<unsigned> nobc_sn_linkage = mesh_no_bc_data->get_side_to_node_linkage();
    FAIL_IF_NOT(std::is_permutation(nobc_sn_linkage.begin(), nobc_sn_linkage.end(),
                                    side_to_node_linkage.begin()));
    FAIL_IF_NOT(mesh_no_bc_data->get_side_set_flag() == nobc_side_set_flag);
  }

  // get the node-to-cell layout generated by the mesh
  const Draco_Mesh::Dual_Layout nc_layout = mesh->get_nc_linkage();

  // make sure the layout size is the number of nodes
  FAIL_IF_NOT(nc_layout.size() == mesh_iface.num_nodes);

  // sanity check each cell vector per node
  for (unsigned node = 0; node < mesh_iface.num_nodes; ++node) {

    // check possible sizes for cell vector per node
    FAIL_IF_NOT(nc_layout.at(node).size() > 0);
    FAIL_IF(nc_layout.at(node).size() == 3);
    FAIL_IF_NOT(nc_layout.at(node).size() <= 4);

    // these check for emergent orthogonal strides in the dual layout
    if (nc_layout.at(node).size() == 4) {

      // check cell stride
      FAIL_IF_NOT(nc_layout.at(node)[1].first == nc_layout.at(node)[0].first + 1);
      FAIL_IF_NOT(nc_layout.at(node)[3].first == nc_layout.at(node)[2].first + 1);
      FAIL_IF_NOT(nc_layout.at(node)[2].first == nc_layout.at(node)[0].first + num_xdir);

      // check node stride
      FAIL_IF_NOT(nc_layout.at(node)[1].second[0] == nc_layout.at(node)[0].second[1]);
      FAIL_IF_NOT(nc_layout.at(node)[1].second[1] == nc_layout.at(node)[3].second[0]);
      FAIL_IF_NOT(nc_layout.at(node)[2].second[0] == nc_layout.at(node)[3].second[1]);
      FAIL_IF_NOT(nc_layout.at(node)[2].second[1] == nc_layout.at(node)[0].second[0]);
    }

    if (nc_layout.at(node).size() == 2) {

      // check cell stride
      const bool has_xdir_stride = nc_layout.at(node)[1].first == nc_layout.at(node)[0].first + 1;
      const bool has_ydir_stride =
          nc_layout.at(node)[1].first == nc_layout.at(node)[0].first + num_xdir;
      FAIL_IF_NOT(has_xdir_stride || has_ydir_stride);

      // check node stride
      const bool has_node_stride =
          nc_layout.at(node)[1].second[1] == nc_layout.at(node)[0].second[0] ||
          nc_layout.at(node)[1].second[0] == nc_layout.at(node)[0].second[1];
      FAIL_IF_NOT(has_node_stride);
    }
  }

  // successful test output
  if (ut.numFails == 0)
    PASSMSG("2D Draco_Mesh tests okay.");
  return;
}

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  rtt_c4::ParallelUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    Insist(rtt_c4::nodes() == 1, "This test only uses 1 PE.");
    spherical_mesh_1d(ut);
    cartesian_mesh_2d(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of mesh/test/tstDraco_Mesh.cc
//------------------------------------------------------------------------------------------------//
