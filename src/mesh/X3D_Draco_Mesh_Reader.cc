//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   mesh/X3D_Draco_Mesh_Reader.cc
 * \author Ryan Wollaeger <wollaeger@lanl.gov>, Kendra Long
 * \date   Thursday, Jul 12, 2018, 08:46 am
 * \brief  X3D_Draco_Mesh_Reader class implementation file.
 * \note   Copyright (C) 2018-2020 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "X3D_Draco_Mesh_Reader.hh"
#include "ds++/DracoStrings.hh"
#include <fstream>
#include <iostream>
#include <set>

namespace rtt_mesh {

//------------------------------------------------------------------------------------------------//
// CONSTRUCTOR
//------------------------------------------------------------------------------------------------//
/*!
 * \brief X3D_Draco_Mesh_Reader constructor.
 *
 * \param[in] filename_ name of file to be parsed
 * \param[in] bdy_filenames_ names of files with lists of side node indexes
 * \param[in] bdy_flags_ unsigned int indicating B.C. per side file (bdy_filenames_)
 */
X3D_Draco_Mesh_Reader::X3D_Draco_Mesh_Reader(const std::string &filename_,
                                             const std::vector<std::string> &bdy_filenames_,
                                             const std::vector<unsigned> &bdy_flags_)
    : filename(filename_), bdy_filenames(bdy_filenames_), bdy_flags(bdy_flags_) {
  // check for valid file name
  Require(filename_.size() > 0);
  Require(bdy_flags_.size() <= bdy_filenames_.size());
}

//------------------------------------------------------------------------------------------------//
// PUBLIC FUNCTIONS
//------------------------------------------------------------------------------------------------//
/*!
 * \brief Build the cell-face index map to the corresponding coordinates.
 */
void X3D_Draco_Mesh_Reader::read_mesh() {

  Require(parsed_pairs.size() == 0);

  // STEP 1: open the file stream

  std::ifstream x3dfile(filename.c_str());

  // file must exist and be readable
  if (!x3dfile.is_open()) {
    throw std::runtime_error("Failed to find or open specified X3D mesh file " + filename);
  }

  // STEP 2: parse file token stream into an initial vector of string pairs

  std::vector<std::pair<std::string, std::string>> raw_pairs;

  while (!x3dfile.eof()) {
    // get a line from the file
    std::string data_line;
    std::getline(x3dfile, data_line);

    // trim outermost white space
    data_line = rtt_dsxx::trim(data_line);

    // ignore empty lines
    if (data_line.size() == 0)
      continue;

    // find last char index of key
    size_t key_offset = data_line.find(' ');

    // first string without whitespace will be key
    std::string key = data_line.substr(0, key_offset);

    // check if anything is after key, to set as value
    std::string value = "";
    if (key_offset != std::string::npos)
      value = data_line.substr(key_offset);

    // add pairing even if value string is size 0 (gets headers and footers)
    raw_pairs.emplace_back(std::pair<std::string, std::string>(key, value));
  }

  // STEP 3: close the file stream

  x3dfile.close();

  // STEP 4: split strings in value at whitespace

  for (auto rpair : raw_pairs) {

    // split the string value into a vector
    std::vector<std::string> tmp_vec = rtt_dsxx::tokenize(rpair.second);

    // add string-vector pair to intermediate map
    parsed_pairs.push_back(Parsed_Element(rpair.first, tmp_vec));
  }

  // STEP 5: derive mesh data maps from parsed_pairs

  // keep track of point parsed_pairs data
  size_t dist = 0;
  if (parsed_pairs[0].first != "header")
    dist++;

  // parse x3d header block and generate x3d_header_map
  Remember(size_t dist_old = dist);
  x3d_header_map = map_x3d_block<std::string, size_t>("header", dist);
  Check(dist > 0);
  Check(dist > dist_old);

  // parse x3d node coordinate block and generate x3d_coord_map
  Remember(dist_old = dist);
  x3d_coord_map = map_x3d_block<int, double>("nodes", dist);
  Check(dist > dist_old);

  // parse x3d face-to-node index map
  Remember(dist_old = dist);
  x3d_facenode_map = map_x3d_block<int, int>("faces", dist);
  Check(dist > dist_old);

  // parse x3d cell-to-face index map
  Remember(dist_old = dist);
  x3d_cellface_map = map_x3d_block<int, int>("cells", dist);
  Check(dist > dist_old);

  // parse x3d material flags
  while (dist < parsed_pairs.size() && parsed_pairs[dist].first != "matid")
    dist++;

  if (parsed_pairs[dist].first != "matid") {
    throw std::invalid_argument("no matid found in x3d file " + filename);
  }

  dist++;
  // As the X3D manual notes, matid is a little weird.
  x3d_matids.resize(get_numcells());
  unsigned i = 0;
  while (i < x3d_matids.size()) {
    if (dist > parsed_pairs.size()) {
      throw std::invalid_argument("wrong number of matids in x3d file " + filename);
    }
    x3d_matids[i++] = parsed_pairs[dist].first;
    auto const &matids = parsed_pairs[dist].second;
    dist++;
    for (auto const &id : matids) {
      x3d_matids[i++] = id;
    }
  }

  // STEP 6: parse side node indices and map to faces

  if (bdy_filenames.size() > 0)
    read_bdy_files();

  Ensure(parsed_pairs.size() > 0);
  Ensure(x3d_header_map.size() > 0);
  Ensure(x3d_coord_map.size() > 0);
  Ensure(x3d_facenode_map.size() > 0);
  Ensure(x3d_cellface_map.size() > 0);
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Return number of nodes for a given cell.
 *
 * \param[in] cell index of cell
 *
 * \return number of nodes for cell
 */
unsigned X3D_Draco_Mesh_Reader::get_celltype(size_t cell) const {

  // x3d file's node, face, and cell indexes start from 1
  Check(cell + 1 < INT_MAX);
  const std::vector<int> &cell_data = x3d_cellface_map.at(static_cast<int>(cell + 1));
  const size_t num_faces = cell_data[0];

  Ensure(num_faces > 0);
  Ensure(num_faces < UINT_MAX);
  return static_cast<unsigned>(num_faces);
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Return the vector of node indices for a given cell.
 *
 * \param[in] cell index of cell
 *
 * \return vector of int node indices
 */
std::vector<unsigned> X3D_Draco_Mesh_Reader::get_cellnodes(size_t cell) const {

  Require(cell < static_cast<size_t>(x3d_header_map.at("elements")[0]));

  // x3d file's node, face, and cell indexes start from 1
  Check(cell + 1 < INT_MAX);
  const std::vector<int> &cell_data = x3d_cellface_map.at(static_cast<int>(cell + 1));
  const size_t num_faces = cell_data[0];

  // calculate number of nodes for this cell
  std::vector<unsigned> node_indexes;

  // track unique node entries with a set
  std::set<int> node_index_set;
  for (size_t i = 1; i <= num_faces; ++i) {

    // get the face index, which will by key for face-to-node map
    int face = cell_data[i];

    // get a vector of nodes for this face
    std::vector<unsigned> tmp_vec = get_facenodes(face);

    // insert into the cell vector
    for (auto j : tmp_vec)
      node_indexes.push_back(j);
  }

  // subtract 1 to get base 0 nodes
  for (auto &node_index : node_indexes)
    node_index--;

  Ensure(node_indexes.size() > 0);
  return node_indexes;
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Return the vector of node indices for a given cell.
 *
 * \param[in] cell index of cell
 * \param[in] face face index local to cell
 *
 * \return vector of int node indices
 */
std::vector<unsigned> X3D_Draco_Mesh_Reader::get_cellfacenodes(size_t cell, size_t face) const {

  Require(cell < static_cast<size_t>(x3d_header_map.at("elements")[0]));

  // x3d file's node, face, and cell indexes start from 1
  Check(cell + 1 < INT_MAX);
  const std::vector<int> &cell_data = x3d_cellface_map.at(static_cast<int>(cell + 1));
  Remember(const size_t num_faces = cell_data[0]);
  Check(face < num_faces);

  // get the face index, which will by key for face-to-node map
  int map_face = cell_data[face + 1];

  // get a vector of nodes for this face
  std::vector<unsigned> node_indexes = get_facenodes(map_face);

  // subtract 1 to get base 0 nodes
  for (auto &node_index : node_indexes)
    node_index--;

  Ensure(node_indexes.size() > 0);
  return node_indexes;
}

//------------------------------------------------------------------------------------------------//
// PRIVATE FUNCTIONS
//------------------------------------------------------------------------------------------------//
/*!
 * \brief Find iterator of a vector of pairs at a key value.
 *
 * \param[in] pairs const reference to the vector of pairs
 * \param[in] key string key to search
 * \param[in] start index from which to start search for key
 *
 * \return iterator to pair with key
 */
X3D_Draco_Mesh_Reader::Parsed_Elements::const_iterator
X3D_Draco_Mesh_Reader::find_iter_of_key(const Parsed_Elements &pairs, std::string key,
                                        size_t start) {
  auto start_it = pairs.begin() + start;
  auto it = std::find_if(start_it, pairs.end(),
                         [&key](const Parsed_Element &p) { return p.first == key; });
  return it;
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Simply return key as specialization to numerical conversion.
 *
 * \param[in] skey string key
 *
 * \return numerical key of type "KT"
 */
template <> std::string X3D_Draco_Mesh_Reader::convert_key<std::string>(const std::string &skey) {
  std::string ret_key = skey;
  return ret_key;
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Return the vector of node indices for a given face.
 *
 * \param[in] face index of face
 *
 * \return vector of int node indices
 */
std::vector<unsigned> X3D_Draco_Mesh_Reader::get_facenodes(size_t face) const {

  Require(face <= x3d_header_map.at("faces")[0]);
  Check(face < INT_MAX);

  // number of nodes is first value after face index in x3d file
  const std::vector<int> &face_data = x3d_facenode_map.at(static_cast<int>(face));
  const size_t num_nodes = face_data[0];

  // return vector
  std::vector<unsigned> node_indexes(num_nodes);

  // push each node instance onto node vector (subtract 1 to get 0-based node)
  for (size_t j = 1; j <= num_nodes; ++j)
    node_indexes[j - 1] = face_data[j];

  Ensure(node_indexes.size() > 0);
  return node_indexes;
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Read side node lists from given boundary files
 */
void X3D_Draco_Mesh_Reader::read_bdy_files() {

  Require(bdy_filenames.size() > 0);
  Require(x3d_header_map.size() > 0);
  Require(x3d_facenode_map.size() > 0);

  const size_t num_flag = bdy_flags.size();
  const size_t num_bdy = bdy_filenames.size();

  size_t bdy_key = 0;
  for (auto bdy_fname : bdy_filenames) {

    // open file
    std::ifstream bdy_file(bdy_fname.c_str());

    // file must exist and be readable
    Insist(bdy_file.is_open(), "Failed to find or open specified X3D mesh file.");

    // append entries to vector of side nodes
    while (!bdy_file.eof()) {

      // get a line from the file
      std::string data_line;
      std::getline(bdy_file, data_line);

      // trim outermost white space
      data_line = rtt_dsxx::trim(data_line);

      // ignore empty lines
      if (data_line.size() == 0)
        continue;

      // try converting to integer
      unsigned side_node;
      try {
        side_node = rtt_dsxx::parse_number_impl<int>(data_line);
      } catch (std::invalid_argument &err) {
        Insist(false, err.what());
      }

      // add to flag-node map
      bc_node_map[bdy_key].push_back(side_node);
    }

    // close the file
    bdy_file.close();

    // increment boundary key
    bdy_key++;
  }

  // Insist that there was at least one side node in all the files
  Insist(bc_node_map.size() > 0, "Boundary file(s) read, but no side nodes.");

  // treat sides as a subset of cell faces here
  int num_side = 0;
  for (size_t bdy = 0; bdy < num_bdy; ++bdy) {

    // calculate flag key and get reference to associated side node vector
    const unsigned flag_key = bdy < num_flag ? bdy_flags[bdy] : 0;
    std::vector<unsigned> &flag_node_vec = bc_node_map[bdy];
    std::sort(flag_node_vec.begin(), flag_node_vec.end());

    // find the mesh faces that have nodes in this flags set
    for (auto const &face_nodes : x3d_facenode_map) {

      // sort vector of nodes associated with this face
      std::vector<unsigned> fnode_vec = get_facenodes(face_nodes.first);
      std::sort(fnode_vec.begin(), fnode_vec.end());

      // \todo: check for node index duplicates

      // find common nodes between side nodes and face
      std::vector<unsigned> nodes_in_common;
      std::set_intersection(flag_node_vec.begin(), flag_node_vec.end(), fnode_vec.begin(),
                            fnode_vec.end(), std::back_inserter(nodes_in_common));

      // if the face is entirely composed of side nodes, then it is a side
      if (nodes_in_common == fnode_vec) {

        // add to the side-node map
        x3d_sidenode_map.insert(
            std::pair<int, std::vector<unsigned>>(num_side, get_facenodes(face_nodes.first)));

        // add to the side-flag map
        x3d_sideflag_map.insert(std::pair<int, unsigned>(num_side, flag_key));

        // increment side counter
        num_side++;
      }
    }
  }

  // decrement node indices
  for (int j = 0; j < num_side; ++j) {
    for (auto &node_index : x3d_sidenode_map.at(j))
      node_index--;
  }
  for (size_t j = 0; j < bc_node_map.size(); ++j) { // NOLINT
    for (auto &node_index : bc_node_map.at(j))
      node_index--;
  }

  Ensure(x3d_sidenode_map.size() > 0);
}

} // end namespace rtt_mesh

//------------------------------------------------------------------------------------------------//
// end of mesh/X3D_Draco_Mesh_Reader.cc
//------------------------------------------------------------------------------------------------//
