//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   mesh/Draco_Mesh_Builder.t.hh
 * \author Ryan Wollaeger <wollaeger@lanl.gov>
 * \date   Tuesday, Jul 03, 2018, 11:26 am
 * \brief  Draco_Mesh_Builder class header file.
 * \note   Copyright (C) 2018-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "Draco_Mesh.hh"
#include "Draco_Mesh_Builder.hh"
#include "ds++/Assert.hh"
#include <algorithm>
#include <iostream>
#include <numeric>

//! \bug this can be fixed when [1] is fixed via C++17 'constexpr if'
#if defined(MSVC)
#pragma warning(push)
#pragma warning(disable : 4702)
#endif

namespace rtt_mesh {

//------------------------------------------------------------------------------------------------//
// CONSTRUCTOR
//------------------------------------------------------------------------------------------------//
/*!
 * \brief Draco_Mesh_Builder constructor.
 *
 * \param[in] reader_ shared pointer to a mesh reader object
 */
template <typename FRT>
Draco_Mesh_Builder<FRT>::Draco_Mesh_Builder(std::shared_ptr<FRT> reader_) : reader(reader_) {
  Require(reader_ != nullptr);
}

//------------------------------------------------------------------------------------------------//
// PUBLIC INTERFACE
//------------------------------------------------------------------------------------------------//
/*!
 * \brief Build a Draco_Mesh object
 *
 * \param[in] geometry enumeration of mesh geometry
 *
 * \return shared pointer to the Draco_Mesh object
 */
template <typename FRT>
std::shared_ptr<Draco_Mesh>
Draco_Mesh_Builder<FRT>::build_mesh(rtt_mesh_element::Geometry geometry) {

  Require(geometry != rtt_mesh_element::Geometry::END_GEOMETRY);

  // >>> GENERATE MESH CONSTRUCTOR ARGUMENTS

  // get the number of dimensions
  unsigned dimension = reader->get_numdim();
  Check(dimension >= 1);
  Check(dimension <= 3);

  // get the number of cells
  size_t num_cells = reader->get_numcells();

  // generate the cell type vector
  std::vector<unsigned> num_faces_per_cell(num_cells);
  for (size_t cell = 0; cell < num_cells; ++cell)
    num_faces_per_cell[cell] = reader->get_celltype(cell);

  // \todo: Can the cell definitions past num_cells - 1 be checked as invalid?

  // generate the cell-to-node linkage
  std::vector<unsigned> cell_to_node_linkage;
  for (size_t cell = 0; cell < num_cells; ++cell) {

    // insert the vector of node indices
    const std::vector<unsigned> cell_nodes = reader->get_cellnodes(cell);
    cell_to_node_linkage.insert(cell_to_node_linkage.end(), cell_nodes.begin(), cell_nodes.end());
  }

  // get the number of sides
  size_t num_sides = reader->get_numsides();

  // generate the side node count vector
  size_t sn_linkage_size = 0;
  std::vector<unsigned> side_node_count(num_sides);
  std::vector<unsigned> side_set_flag(num_sides);
  for (size_t side = 0; side < num_sides; ++side) {

    // acquire the number of nodes associated with this side def
    Check(reader->get_sidetype(side) < UINT_MAX);
    side_node_count[side] = static_cast<unsigned>(reader->get_sidetype(side));

    // this is not required in RTT meshes, but is so in Draco_Mesh
    Check(dimension == 2 ? side_node_count[side] == 2 : true);

    // get the 1st side flag associated with this side
    // \todo: What happens when side has no flags?
    side_set_flag[side] = reader->get_sideflag(side);

    // increment size of cell-to-node linkage array
    sn_linkage_size += side_node_count[side];
  }

  // \todo: Can the side definitions past num_sides - 1 be checked as invalid?

  // generate the side-to-node linkage
  std::vector<unsigned> side_to_node_linkage;
  side_to_node_linkage.reserve(sn_linkage_size);
  for (size_t side = 0; side < num_sides; ++side) {

    // insert the vector of node indices
    const std::vector<unsigned> side_nodes = reader->get_sidenodes(side);
    side_to_node_linkage.insert(side_to_node_linkage.end(), side_nodes.begin(), side_nodes.end());
  }

  // get the number of nodes
  size_t num_nodes = reader->get_numnodes();

  Check(num_nodes >= num_cells);

  // generate the global node number serialized vector of coordinates
  std::vector<unsigned> global_node_number(num_nodes);
  std::vector<double> coordinates(dimension * num_nodes);
  for (size_t node = 0; node < num_nodes; ++node) {

    // set the "global" node indices
    Check(node < UINT_MAX);
    global_node_number[node] = static_cast<unsigned>(node);

    // get coordinates for this node
    const std::vector<double> node_coord = reader->get_nodecoord(node);

    // populate coordinate vector
    for (unsigned d = 0; d < dimension; ++d)
      coordinates[dimension * node + d] = node_coord[d];
  }

  // reserve some space for num_nodes_per_face_per_cell
  std::vector<unsigned> num_nodes_per_face_per_cell;
  num_nodes_per_face_per_cell.reserve(
      std::accumulate(num_faces_per_cell.begin(), num_faces_per_cell.end(), 0U));

  // generate num_nodes_per_face_per_cell vector
  unsigned cf_counter = 0;
  for (size_t cell = 0; cell < num_cells; ++cell) {
    for (unsigned face = 0; face < num_faces_per_cell[cell]; ++face) {

      // store number of nodes for this face
      num_nodes_per_face_per_cell.push_back(
          static_cast<unsigned>(reader->get_cellfacenodes(cell, face).size()));

      // increment counter
      cf_counter++;
    }
  }

  Remember(auto cn_minmax =
               std::minmax_element(cell_to_node_linkage.begin(), cell_to_node_linkage.end()));
  Remember(auto sn_minmax =
               std::minmax_element(side_to_node_linkage.begin(), side_to_node_linkage.end()));
  Ensure(*cn_minmax.second < num_nodes);
  Ensure(side_to_node_linkage.size() > 0 ? *sn_minmax.second < num_nodes : true);

  // >>> CONSTRUCT THE MESH

  std::shared_ptr<Draco_Mesh> mesh(new Draco_Mesh(
      dimension, geometry, num_faces_per_cell, cell_to_node_linkage, side_set_flag, side_node_count,
      side_to_node_linkage, coordinates, global_node_number, num_nodes_per_face_per_cell));

  return mesh;
}

} // end namespace rtt_mesh

#if defined(MSVC)
#pragma warning(pop)
#endif

//------------------------------------------------------------------------------------------------//
// end of mesh/Draco_Mesh_Builder.t.hh
//------------------------------------------------------------------------------------------------//
