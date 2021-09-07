//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   mesh/Draco_Mesh.hh
 * \author Ryan Wollaeger <wollaeger@lanl.gov>
 * \date   Thursday, Jun 07, 2018, 15:38 pm
 * \brief  Draco_Mesh class header file.
 * \note   Copyright (C) 2018-2021 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_mesh_Draco_Mesh_hh
#define rtt_mesh_Draco_Mesh_hh

#include "ds++/config.h"
#include "mesh_element/Geometry.hh"
#include <array>
#include <map>
#include <set>
#include <vector>

namespace rtt_mesh {

//================================================================================================//
/*!
 * \class Draco_Mesh
 *
 * \brief General unstructured mesh class.
 *
 * The Draco_Mesh class takes cell-node (or cell-vertex) data, and generates a mesh with layout
 * (cell adjacency) information.  This class also provides basic services, including access to cell
 * information.  This mesh is based on an unstructured mesh implementation by Kent Budge.
 *
 * Important features for a fully realized Draco_Mesh are the following:
 * 1) Geometry, which implies a metric for distance between points.
 * 2) Layout, which stores cell connectivity and hence the mesh topology.
 *    a) It has an internal layout containing local cell-to-cell linkage,
 *    b) a boundary layout with true-boundary linkage, and
 *    c) a ghost layout containing cell-to-ghost-cell linkage.
 * 3) Dual_Layout, which stores node connectivity and effectively inverts Layout
 * 4) Dual_Ghost_Layout, which stores node connectivity to off-process adjacent cells and nodes.
 *    This an has additional field for the MPI rank index the neighboring cell and nodes are on.
 *
 * Possibly temporary features:
 * 1) The num_faces_per_cell_ vector (argument to the constructor) is currently taken to be the
 *    number of faces per cell.
 * 2) The layout data structure(s) will probably be moved to a separate class, where accessors might
 *    be used on a flattened version.
 */
//================================================================================================//

class Draco_Mesh {
public:
  // >>> TYPEDEFS

  using Geometry = rtt_mesh_element::Geometry;

  // e.g.: (key: cell, value: vector of pairs of neighbors and shared nodes per face)
  using Layout =
      std::map<unsigned int, std::vector<std::pair<unsigned int, std::vector<unsigned int>>>>;

  // e.g.: (key: node, value: vector of local or ghost cell indices), only 2D
  using CellNodes_Pair = std::pair<unsigned int, std::array<unsigned int, 2>>;
  using Dual_Layout = std::map<unsigned int, std::vector<CellNodes_Pair>>;

  // e.g.: (key: node, value: vector of pairs of rank and local cell index on the rank)
  // vectors of pairs are ordered in increasing rank, by construction (see Draco_Mesh.cc)
  using Coord_NBRS = std::pair<std::array<double, 2>, std::array<double, 2>>;
  using Dual_Ghost_Layout =
      std::map<unsigned int, std::vector<std::pair<CellNodes_Pair, unsigned int>>>;
  using Dual_Ghost_Layout_Coords = std::map<unsigned int, std::vector<Coord_NBRS>>;

protected:
  // >>> DATA

  // Dimension
  const unsigned dimension;

  // Geometry enumeration
  const Geometry geometry;

  // Number of cells
  const unsigned num_cells;

  // Number of nodes
  const unsigned num_nodes;

  // Side set flag (can be used for mapping BCs to sides)
  std::vector<unsigned> side_set_flag;

  // Ghost cell indices local to a different node, subscripted with a local ghost cell index
  const std::vector<int> ghost_cell_number;

  // Node index for each ghost cell, subscripted with local ghost cell index
  const std::vector<int> ghost_cell_rank;

  // Vector subscripted with node index with coordinate vector
  const std::vector<std::vector<double>> node_coord_vec;

  // Cell types and node indices per cell
  const std::vector<unsigned> m_num_faces_per_cell;
  const std::vector<unsigned> m_num_nodes_per_face_per_cell;
  const std::vector<std::vector<std::vector<unsigned>>> m_cell_to_node_linkage;

  // Side types and node indices per side
  std::vector<unsigned> m_side_node_count;
  std::vector<unsigned> m_side_to_node_linkage;

  // Layout of mesh: vector index is cell index, vector element is description of cell's adjacency
  // to other cells in the mesh.
  Layout cell_to_cell_linkage;

  // Side layout of mesh
  Layout cell_to_side_linkage;

  // Ghost cell layout of mesh
  Layout cell_to_ghost_cell_linkage;

  // Node map to vector of local cells
  Dual_Layout node_to_cellnode_linkage;

  // Node map to vector of ghost cells
  Dual_Ghost_Layout node_to_ghost_cell_linkage;

  // Node map to vector of adjacent coordinates bounding adjacent ghost cells
  Dual_Ghost_Layout_Coords node_to_ghost_coord_linkage;

public:
  //! Constructor.
  Draco_Mesh(unsigned dimension_, Geometry geometry_,
             const std::vector<unsigned> &num_faces_per_cell_,
             const std::vector<unsigned> &cell_to_node_linkage_,
             const std::vector<unsigned> side_set_flag_,
             const std::vector<unsigned> &side_node_count_,
             const std::vector<unsigned> &side_to_node_linkage_,
             const std::vector<double> &coordinates_,
             const std::vector<unsigned> &global_node_number_,
             const std::vector<unsigned> &num_nodes_per_face_per_cell_,
             const std::vector<unsigned> &ghost_cell_type_ = {},
             const std::vector<unsigned> &ghost_cell_to_node_linkage_ = {},
             const std::vector<int> &ghost_cell_number_ = {},
             const std::vector<int> &ghost_cell_rank_ = {});

  // >>> ACCESSORS

  unsigned get_dimension() const { return dimension; }
  Geometry get_geometry() const { return geometry; }
  unsigned get_num_cells() const { return num_cells; }
  unsigned get_num_nodes() const { return num_nodes; }
  std::vector<unsigned> get_side_set_flag() const { return side_set_flag; }
  std::vector<int> get_ghost_cell_numbers() const { return ghost_cell_number; }
  std::vector<int> get_ghost_cell_ranks() const { return ghost_cell_rank; }
  std::vector<std::vector<double>> get_node_coord_vec() const { return node_coord_vec; }
  std::vector<unsigned> get_num_faces_per_cell() const { return m_num_faces_per_cell; }
  std::vector<unsigned> get_num_nodes_per_face_per_cell() const {
    return m_num_nodes_per_face_per_cell;
  }
  std::vector<std::vector<std::vector<unsigned>>> get_cell_to_node_linkage() const {
    return m_cell_to_node_linkage;
  }
  std::vector<unsigned> get_side_node_count() const { return m_side_node_count; }
  std::vector<unsigned> get_side_to_node_linkage() const { return m_side_to_node_linkage; }
  Layout get_cc_linkage() const { return cell_to_cell_linkage; }
  Layout get_cs_linkage() const { return cell_to_side_linkage; }
  Layout get_cg_linkage() const { return cell_to_ghost_cell_linkage; }
  Dual_Layout get_nc_linkage() const { return node_to_cellnode_linkage; }
  Dual_Ghost_Layout get_ngc_linkage() const { return node_to_ghost_cell_linkage; }
  Dual_Ghost_Layout_Coords get_ngcoord_linkage() const { return node_to_ghost_coord_linkage; }

  // >>> SERVICES

  const std::vector<unsigned> get_cell_nodes(const unsigned cell) const;
  const std::vector<unsigned> get_flat_cell_node_linkage() const;

private:
  // >>> SUPPORT FUNCTIONS

  //! Calculate (merely deserialize) the vector of node coordinates
  std::vector<std::vector<double>>
  compute_node_coord_vec(const std::vector<double> &coordinates) const;

  //! Convert cell-node linkage into tensor for easier indexing by cell, face
  std::vector<std::vector<std::vector<unsigned>>>
  compute_cell_to_node_tensor(const std::vector<unsigned> &num_faces_per_cell,
                              const std::vector<unsigned> &num_nodes_per_face_per_cell,
                              const std::vector<unsigned> &cell_to_node_linkage) const;

  //! Calculate the cell-to-cell linkage with face type vector
  void compute_cell_to_cell_linkage(const std::vector<unsigned> &num_faces_per_cell,
                                    const std::vector<unsigned> &cell_to_node_linkage,
                                    const std::vector<unsigned> &num_nodes_per_face_per_cell,
                                    const std::vector<unsigned> &side_node_count,
                                    const std::vector<unsigned> &side_to_node_linkage,
                                    const std::vector<unsigned> &ghost_cell_type,
                                    const std::vector<unsigned> &ghost_cell_to_node_linkage);

  //! Calculate a map of node vectors to indices (sides, ghost cells)
  std::map<std::set<unsigned>, unsigned>
  compute_node_vec_indx_map(const std::vector<unsigned> &indx_type,
                            const std::vector<unsigned> &indx_to_node_linkage) const;

  //! Calculate cell-corner-cell layouts (adjacent cells not sharing a face)
  void compute_node_to_cell_linkage(const std::vector<unsigned> &ghost_cell_type,
                                    const std::vector<unsigned> &ghost_cell_to_node_linkage,
                                    const std::vector<unsigned> &global_node_number);
};

} // end namespace rtt_mesh

#endif // rtt_mesh_Draco_Mesh_hh

//------------------------------------------------------------------------------------------------//
// end of mesh/Draco_Mesh.hh
//------------------------------------------------------------------------------------------------//
