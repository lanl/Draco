//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   RTT_Format_Reader/Dims.cc
 * \author B.T. Adams
 * \date   Wed Jun 7 10:33:26 2000
 * \brief  Implementation file for RTT_Format_Reader/Dims class.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "Dims.hh"
#include "ds++/Assert.hh"

namespace rtt_RTT_Format_Reader {
/*!
 * \brief Parses the dims (dimensions) data block from the mesh file via calls to private member
 *        functions.
 * \param meshfile Mesh file name.
 */
void Dims::readDims(ifstream &meshfile) {
  readKeyword(meshfile);
  readUnits(meshfile);
  readCellDefs(meshfile);
  readDimensions(meshfile);
  readNodes(meshfile);
  readSides(meshfile);
  readCells(meshfile);
  readEndKeyword(meshfile);
}
/*!
 * \brief Reads and validates the dims (dimensions) block keyword.
 * \param meshfile Mesh file name.
 */
void Dims::readKeyword(ifstream &meshfile) {
  string dummyString;

  meshfile >> dummyString;
  Insist(dummyString == "dims", "Invalid mesh file: Dimension block missing");
  std::getline(meshfile, dummyString); // read and discard blank line.
}
/*!
 * \brief Reads and validates the dims (dimensions) coordinates and time units.
 * \param meshfile Mesh file name.
 */
void Dims::readUnits(ifstream &meshfile) {
  string dummyString;

  meshfile >> dummyString >> coor_units;
  Insist(dummyString == "coor_units", "Invalid mesh file: Dimension block missing coor_units");
  std::getline(meshfile, dummyString);

  meshfile >> dummyString >> prob_time_units;
  Insist(dummyString == "prob_time_units",
         "Invalid mesh file: Dimension block missing prob_time_units");
  std::getline(meshfile, dummyString);
}
/*!
 * \brief Reads and validates the dims (dimensions) cell definition data.
 * \param meshfile Mesh file name.
 */
void Dims::readCellDefs(ifstream &meshfile) {
  string dummyString;

  meshfile >> dummyString >> ncell_defs;
  Insist(dummyString == "ncell_defs", "Invalid mesh file: Dimension block missing ncell_defs");
  std::getline(meshfile, dummyString);

  meshfile >> dummyString >> nnodes_max;
  Insist(dummyString == "nnodes_max", "Invalid mesh file: Dimension block missing nnodes_max");
  std::getline(meshfile, dummyString);

  meshfile >> dummyString >> nsides_max;
  Insist(dummyString == "nsides_max", "Invalid mesh file: Dimension block missing nsides_max");
  std::getline(meshfile, dummyString);

  meshfile >> dummyString >> nnodes_side_max;
  Insist(dummyString == "nnodes_side_max",
         "Invalid mesh file: Dimension block missing nnodes_side_max");
  std::getline(meshfile, dummyString);
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Reads and validates the dims (dimensions) spatial and topological dimension data.
 * \param meshfile Mesh file name.
 */
void Dims::readDimensions(ifstream &meshfile) {
  string dummyString;

  meshfile >> dummyString >> ndim;
  Insist(dummyString == "ndim", "Invalid mesh file: Dimension block missing ndim");
  std::getline(meshfile, dummyString);

  meshfile >> dummyString >> ndim_topo;
  Insist(dummyString == "ndim_topo", "Invalid mesh file: Dimension block missing ndim_topo");
  std::getline(meshfile, dummyString);
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Reads and validates the dims (dimensions) node data.
 * \param meshfile Mesh file name.
 */
void Dims::readNodes(ifstream &meshfile) {
  string dummyString;

  meshfile >> dummyString >> nnodes;
  Insist(dummyString == "nnodes", "Invalid mesh file: Dimension block missing nnodes");
  std::getline(meshfile, dummyString);

  meshfile >> dummyString >> nnode_flag_types;
  Insist(dummyString == "nnode_flag_types",
         "Invalid mesh file: Dimension block missing nnode_flag_types");
  std::getline(meshfile, dummyString);

  nnode_flags.resize(nnode_flag_types);
  meshfile >> dummyString;
  Insist(dummyString == "nnode_flags", "Invalid mesh file: Dimension block missing nnode_flags");
  for (auto iter = nnode_flags.begin(); iter < nnode_flags.end(); ++iter)
    meshfile >> *iter;
  std::getline(meshfile, dummyString);

  meshfile >> dummyString >> nnode_data;
  Insist(dummyString == "nnode_data", "Invalid mesh file: Dimension block missing nnode_data");
  std::getline(meshfile, dummyString);
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Reads and validates the dims (dimensions) side data.
 * \param meshfile Mesh file name.
 */
void Dims::readSides(ifstream &meshfile) {
  string dummyString;

  meshfile >> dummyString >> nsides;
  Insist(dummyString == "nsides", "Invalid mesh file: Dimension block missing nsides");
  std::getline(meshfile, dummyString);

  meshfile >> dummyString >> nside_types;
  Insist(dummyString == "nside_types", "Invalid mesh file: Dimension block missing nside_types");
  std::getline(meshfile, dummyString);

  side_types.resize(nside_types);
  meshfile >> dummyString;
  Insist(dummyString == "side_types", "Invalid mesh file: Dimension block missing side_types");
  for (auto iter = side_types.begin(); iter < side_types.end(); ++iter) {
    meshfile >> *iter;
    --(*iter);
  }
  std::getline(meshfile, dummyString);

  meshfile >> dummyString >> nside_flag_types;
  Insist(dummyString == "nside_flag_types",
         "Invalid mesh file: Dimension block missing nside_flag_types");
  std::getline(meshfile, dummyString);

  nside_flags.resize(nside_flag_types);
  meshfile >> dummyString;
  Insist(dummyString == "nside_flags", "Invalid mesh file: Dimension block missing nside_flags");
  for (auto iter = nside_flags.begin(); iter < nside_flags.end(); ++iter)
    meshfile >> *iter;
  std::getline(meshfile, dummyString);

  meshfile >> dummyString >> nside_data;
  Insist(dummyString == "nside_data", "Invalid mesh file: Dimension block missing nside_data");
  std::getline(meshfile, dummyString);
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Reads and validates the dims (dimensions) cell data.
 * \param meshfile Mesh file name.
 */
void Dims::readCells(ifstream &meshfile) {
  string dummyString;

  meshfile >> dummyString >> ncells;
  Insist(dummyString == "ncells", "Invalid mesh file: Dimension block missing ncells");
  std::getline(meshfile, dummyString);

  meshfile >> dummyString >> ncell_types;
  Insist(dummyString == "ncell_types", "Invalid mesh file: Dimension block missing ncell_types");
  std::getline(meshfile, dummyString);

  cell_types.resize(ncell_types);
  meshfile >> dummyString;
  Insist(dummyString == "cell_types", "Invalid mesh file: Dimension block missing cell_types");
  for (auto iter = cell_types.begin(); iter < cell_types.end(); ++iter) {
    meshfile >> *iter;
    --(*iter);
  }
  std::getline(meshfile, dummyString);

  meshfile >> dummyString >> ncell_flag_types;
  Insist(dummyString == "ncell_flag_types",
         "Invalid mesh file: Dimension block missing ncell_flag_types");
  std::getline(meshfile, dummyString);

  ncell_flags.resize(ncell_flag_types);
  meshfile >> dummyString;
  Insist(dummyString == "ncell_flags", "Invalid mesh file: Dimension block missing ncell_flags");
  for (auto iter = ncell_flags.begin(); iter < ncell_flags.end(); ++iter)
    meshfile >> *iter;
  std::getline(meshfile, dummyString);

  meshfile >> dummyString >> ncell_data;
  Insist(dummyString == "ncell_data", "Invalid mesh file: Dimension block missing ncell_data");
  std::getline(meshfile, dummyString);
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Reads and validates the end_dims block keyword.
 * \param meshfile Mesh file name.
 */
void Dims::readEndKeyword(ifstream &meshfile) {
  string dummyString;

  meshfile >> dummyString;
  Insist(dummyString == "end_dims", "Invalid mesh file: Dimension block missing end_dims");
  std::getline(meshfile, dummyString); // read and discard blank line.
}

} // end namespace rtt_RTT_Format_Reader

//------------------------------------------------------------------------------------------------//
// end of RTT_Format_Reader/end of Dims.cc
//------------------------------------------------------------------------------------------------//
