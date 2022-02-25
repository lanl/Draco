//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   mesh/RTT_Draco_Mesh_Reader.cc
 * \author Ryan Wollaeger <wollaeger@lanl.gov>
 * \date   Friday, Jul 13, 2018, 08:38 am
 * \brief  RTT_Draco_Mesh_Reader header file.
 * \note   Copyright (C) 2018-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "RTT_Draco_Mesh_Reader.hh"

namespace rtt_mesh {

//------------------------------------------------------------------------------------------------//
// CONSTRUCTOR
//------------------------------------------------------------------------------------------------//
/*!
 * \brief RTT_Draco_Mesh_Reader constructor.
 *
 * \param[in] filename_ name of file to be parsed
 */
RTT_Draco_Mesh_Reader::RTT_Draco_Mesh_Reader(const std::string &filename_) : filename(filename_) {
  // check for valid file name
  Insist(filename_.size() > 0, "No file name supplied.");

  // \todo: remove read_mesh function and read in constructor(?)
  Require(rtt_reader == nullptr);
}

//------------------------------------------------------------------------------------------------//
// PUBLIC FUNCTIONS
//------------------------------------------------------------------------------------------------//
/*!
 * \brief Read the mesh by constructing an RTT_Format_Reader
 */
void RTT_Draco_Mesh_Reader::read_mesh() {

  // insist the file can be opened before trying read_mesh() wrapper
  std::ifstream rttfile(filename.c_str());
  Insist(rttfile.is_open(), "Failed to find or open specified RTT mesh file.");
  rttfile.close();

  rtt_reader = std::make_shared<rtt_RTT_Format_Reader::RTT_Format_Reader>(filename);
  Insist(rtt_reader->get_dims_ndim() < 3, "RTT reader mesh dim >= 3.");
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Get the number of nodes for a cell.
 *
 * \param[in] cell index of cell (0-based)
 *
 * \return number of nodes for cell
 */
unsigned RTT_Draco_Mesh_Reader::get_celltype(size_t cell) const {

  // first obtain a cell definition index
  size_t cell_def = rtt_reader->get_cells_type(cell);

  // for Draco_Mesh, cell_type is number of nodes (in 1-2D, this is # of faces)
  Check(rtt_reader->get_cell_defs_nnodes(cell_def) < UINT_MAX);
  auto cell_type = static_cast<unsigned>(rtt_reader->get_cell_defs_nnodes(cell_def));

  return cell_type;
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Get the vector of node indices for a side.
 *
 * \param[in] side index of side (0-based)
 *
 * \return number of nodes for side
 */
size_t RTT_Draco_Mesh_Reader::get_sidetype(size_t side) const {

  // first obtain a side definition index
  size_t side_def = rtt_reader->get_sides_type(side);

  // acquire the number of nodes associated with this side def
  size_t side_type = rtt_reader->get_cell_defs_nnodes(side_def);

  return side_type;
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Get the vector of node indices for a cell.
 *
 * \param[in] cell index of cell (0-based)
 *
 * \return vector of the cell's nodes per face (hence with duplicates)
 */
std::vector<unsigned> RTT_Draco_Mesh_Reader::get_cellnodes(size_t cell) const {

  // get the number of dimensions
  unsigned num_dim = rtt_reader->get_dims_ndim();

  // get the node list for this cell
  const std::vector<unsigned> &cell_node = rtt_reader->get_cells_nodes(cell);

  // get the number of nodes
  size_t cell_type = cell_node.size();

  // internal nodes are not supported in 1D
  Check(num_dim == 1 ? cell_type == 2 : cell_type > 2);

  // vector of node per (cell,face)
  std::vector<unsigned> cellface_node(num_dim * cell_type);

  // set node upper limit based on dimension
  size_t node_limit = num_dim == 1 ? cell_type - 1 : cell_type;
  for (size_t i = 0; i < node_limit; ++i) {

    // get next face index
    size_t j = (i + 1) % cell_type;

    // populate node indices per face
    cellface_node[num_dim * i] = cell_node[i];
    cellface_node[num_dim * i + 1] = cell_node[j];
  }

  return cellface_node;
}

//------------------------------------------------------------------------------------------------//
std::vector<unsigned> RTT_Draco_Mesh_Reader::get_cellfacenodes(size_t cell, size_t face) const {

  // get the number of dimensions
  unsigned const num_dim = rtt_reader->get_dims_ndim();

  // get the node list for this cell
  const std::vector<unsigned> &cell_node = rtt_reader->get_cells_nodes(cell);

  // get the number of nodes
  size_t const cell_type = cell_node.size();

  // internal nodes are not supported in 1D
  Check(num_dim == 1 ? cell_type == 2 : cell_type > 2);

  std::vector<unsigned> face_node;
  face_node.emplace_back(cell_node[face]);
  if (num_dim == 2)
    face_node.emplace_back(cell_node[(face + 1) % cell_type]);

  return face_node;
}

} // end namespace rtt_mesh

//------------------------------------------------------------------------------------------------//
// end of mesh/RTT_Draco_Mesh_Reader.cc
//------------------------------------------------------------------------------------------------//
