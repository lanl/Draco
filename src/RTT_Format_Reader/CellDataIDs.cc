//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   RTT_Format_Reader/CellDataIDs.cc
 * \author B.T. Adams
 * \date   Wed Jun 7 10:33:26 2000
 * \brief  Implementation file for RTT_Format_Reader/CellDataIDs class.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "CellDataIDs.hh"

namespace rtt_RTT_Format_Reader {
/*!
 * \brief Parses the cell_data_ids data block from the mesh file via calls to private member
 *        functions.
 * \param meshfile Mesh file name.
 */
void CellDataIDs::readDataIDs(ifstream &meshfile) {
  readKeyword(meshfile);
  readData(meshfile);
  readEndKeyword(meshfile);
}
/*!
 * \brief Reads and validates the cell_data_ids block keyword.
 * \param meshfile Mesh file name.
 */
void CellDataIDs::readKeyword(ifstream &meshfile) {
  string dummyString;

  meshfile >> dummyString;
  Insist(dummyString == "cell_data_ids", "Invalid mesh file: cell_data_ids block missing");
  std::getline(meshfile, dummyString);
}
/*!
 * \brief Reads and validates the cell_data_ids block data.
 * \param meshfile Mesh file name.
 */
void CellDataIDs::readData(ifstream &meshfile) {
  int dataIDNum;
  string dummyString;

  for (size_t i = 0; i < dims.get_ncell_data(); ++i) {
    Check(i < names.size() && i < units.size());
    meshfile >> dataIDNum >> names[i] >> units[i];
    Insist(static_cast<size_t>(dataIDNum) == i + 1, "Invalid mesh file: cell data ID out of order");
    std::getline(meshfile, dummyString);
  }
}
/*!
 * \brief Reads and validates the end_cell_data_ids block keyword.
 * \param meshfile Mesh file name.
 */
void CellDataIDs::readEndKeyword(ifstream &meshfile) {
  string dummyString;

  meshfile >> dummyString;
  Insist(dummyString == "end_cell_data_ids", "Invalid mesh file: cell_data_ids block missing end");
  std::getline(meshfile, dummyString); // read and discard blank line.
}

} // end namespace rtt_RTT_Format_Reader

//------------------------------------------------------------------------------------------------//
// end of RTT_Format_Reader/CellDataIDs.cc
//------------------------------------------------------------------------------------------------//
