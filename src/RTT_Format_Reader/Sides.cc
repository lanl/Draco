//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   RTT_Format_Reader/Sides.cc
 * \author B.T. Adams
 * \date   Wed Jun 7 10:33:26 2000
 * \brief  Implementation file for RTT_Format_Reader/Sides class.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "Sides.hh"

namespace rtt_RTT_Format_Reader {

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Parses the sides block data from the mesh file via calls to private member functions.
 * \param meshfile Mesh file name.
 */
void Sides::readSides(ifstream &meshfile) {
  readKeyword(meshfile);
  readData(meshfile);
  readEndKeyword(meshfile);
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Reads and validates the sides block keyword.
 * \param meshfile Mesh file name.
 */
void Sides::readKeyword(ifstream &meshfile) {
  string dummyString;

  meshfile >> dummyString;
  Insist(dummyString == "sides", "Invalid mesh file: sides block missing");
  std::getline(meshfile, dummyString);
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Reads and validates the sides block data.
 * \param meshfile Mesh file name.
 */
void Sides::readData(ifstream &meshfile) {
  string dummyString;
  int sideNum;

  for (unsigned i = 0; i < static_cast<unsigned int>(dims.get_nsides()); ++i) {
    sideNum = Nodes::readNextInt(meshfile);
    // meshfile >> sideNum;
    Insist(static_cast<unsigned int>(sideNum) == i + 1,
           "Invalid mesh file: side index out of order");
    Check(i < sideType.size());
    meshfile >> sideType[i];
    --sideType[i];
    Insist(dims.allowed_side_type(sideType[i]), "Invalid mesh file: illegal side type");
    Check(i < nodes.size());
    nodes[i].resize(cellDefs.get_nnodes(sideType[i]));
    for (unsigned j = 0; j < cellDefs.get_nnodes(sideType[i]); ++j) {
      Check(j < nodes[i].size());
      meshfile >> nodes[i][j];
      --nodes[i][j];
    }
    for (unsigned j = 0; j < static_cast<unsigned int>(dims.get_nside_flag_types()); ++j) {
      Check(j < flags[i].size());
      meshfile >> flags[i][j];
      Insist(sideFlags.allowed_flag(j, flags[i][j]), "Invalid mesh file: illegal side flag");
    }
    std::getline(meshfile, dummyString);
  }
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Reads and validates the end_sides block keyword.
 * \param meshfile Mesh file name.
 */
void Sides::readEndKeyword(ifstream &meshfile) {
  string dummyString;

  meshfile >> dummyString;
  Insist(dummyString == "end_sides", "Invalid mesh file: sides block missing end");
  std::getline(meshfile, dummyString); // read and discard blank line.
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Changes the side nodes when the cell definitions specified in the RTT_Format file have
 *        been transformed into an alternative cell definition (e.g., CYGNUS).
 */
void Sides::redefineSides() {
  vector_uint temp_nodes;
  for (size_t st = 0; st < dims.get_nside_types(); st++) {
    int this_side_type = dims.get_side_types(st);
    vector_uint node_map(cellDefs.get_node_map(this_side_type));
    Insist(node_map.size() == cellDefs.get_nnodes(this_side_type), "Error in Sides redefinition.");
    // Check to see if the nodes need to be rearranged for this side type.
    bool redefined = false;
    for (unsigned n = 0; n < node_map.size(); n++) {
      if (node_map[n] != n)
        redefined = true;
    }
    if (redefined) {
      temp_nodes.resize(cellDefs.get_nnodes(this_side_type));
      for (size_t s = 0; s < dims.get_nsides(); s++) {
        if (sideType[s] == this_side_type) {
          for (size_t n = 0; n < nodes[s].size(); n++)
            temp_nodes[node_map[n]] = nodes[s][n];
          nodes[s] = temp_nodes;
        }
      }
    }
    node_map.resize(0);
  }
}

} // end namespace rtt_RTT_Format_Reader

//------------------------------------------------------------------------------------------------//
// end of RTT_Format_Reader/Sides.cc
//------------------------------------------------------------------------------------------------//
