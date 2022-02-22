//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   RTT_Format_Reader/Cells.cc
 * \author B.T. Adams
 * \date   Wed Jun 7 10:33:26 2000
 * \brief  Implementation file for RTT_Format_Reader/Cells class.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "Cells.hh"

namespace rtt_RTT_Format_Reader {

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Parses the cells block data from the mesh file via calls to private member functions.
  * \param meshfile Mesh file name.
 */
void Cells::readCells(ifstream &meshfile) {
  readKeyword(meshfile);
  readData(meshfile);
  readEndKeyword(meshfile);
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Reads and validates the cells block keyword.
 * \param meshfile Mesh file name.
 */
void Cells::readKeyword(ifstream &meshfile) {
  string dummyString;

  meshfile >> dummyString;
  Insist(dummyString == "cells", "Invalid mesh file: cells block missing");
  std::getline(meshfile, dummyString);
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Reads and validates the cells block data.
 * \param meshfile Mesh file name.
 */
void Cells::readData(ifstream &meshfile) {
  string dummyString;
  int cellNum;

  for (size_t i = 0; i < dims.get_ncells(); ++i) {
    cellNum = Nodes::readNextInt(meshfile);
    // meshfile >> cellNum;
    Insist(static_cast<size_t>(cellNum) == i + 1, "Invalid mesh file: cell index out of order");
    Check(i < cellType.size());
    meshfile >> cellType[i];
    --cellType[i];
    Insist(dims.allowed_cell_type(cellType[i]), "Invalid mesh file: illegal cell type");
    Check(i < nodes.size());
    nodes[i].resize(cellDefs.get_nnodes(cellType[i]));
    for (unsigned j = 0; j < cellDefs.get_nnodes(cellType[i]); ++j) {
      Check(j < nodes[i].size());
      meshfile >> nodes[i][j];
      --nodes[i][j];
    }

    for (size_t j = 0; j < dims.get_ncell_flag_types(); ++j) {
      Check(j < flags[i].size());
      meshfile >> flags[i][j];
      Check(j < INT_MAX);
      Insist(cellFlags.allowed_flag(static_cast<int>(j), flags[i][j]),
             "Invalid mesh file: illegal cell flag");
    }
    std::getline(meshfile, dummyString);
  }
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Reads and validates the end_cells block keyword.
 * \param meshfile Mesh file name.
 */
void Cells::readEndKeyword(ifstream &meshfile) {
  string dummyString;

  meshfile >> dummyString;
  Insist(dummyString == "end_cells", "Invalid mesh file: cells block missing end");
  std::getline(meshfile, dummyString); // read and discard blank line.
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Changes the cell nodes when the cell definitions specified in the RTT_Format file have
 *        been transformed into an alternative cell definition (e.g., CYGNUS).
 */
void Cells::redefineCells() {
  vector_uint temp_nodes;
  for (size_t ct = 0; ct < dims.get_ncell_types(); ct++) {
    int this_cell_type = dims.get_cell_types(ct);
    vector_uint node_map(cellDefs.get_node_map(this_cell_type));
    Insist(node_map.size() == cellDefs.get_nnodes(this_cell_type), "Error in Cells redefinition.");

    // Check to see if the nodes need to be rearranged for this cell type.
    bool redefined = false;
    for (size_t n = 0; n < node_map.size(); n++) {
      if (static_cast<size_t>(node_map[n]) != n)
        redefined = true;
    }

    if (redefined) {
      temp_nodes.resize(cellDefs.get_nnodes(this_cell_type));
      for (size_t c = 0; c < dims.get_ncells(); c++) {
        if (cellType[c] == this_cell_type) {
          for (size_t n = 0; n < nodes[c].size(); n++)
            temp_nodes[node_map[n]] = nodes[c][n];
          nodes[c] = temp_nodes;
        }
      }
    }
    node_map.resize(0);
  }
}

} // end namespace rtt_RTT_Format_Reader

//------------------------------------------------------------------------------------------------//
// end of RTT_Format_Reader/Cells.cc
//------------------------------------------------------------------------------------------------//
