//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   RTT_Format_Reader/CellData.cc
 * \author B.T. Adams
 * \date   Wed Jun 7 10:33:26 2000
 * \brief  Implementation file for RTT_Format_Reader/CellData class.
 * \note   Copyright (C) 2010-2021 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "CellData.hh"

namespace rtt_RTT_Format_Reader {

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Parses the cell_data block data from the mesh file via calls to private member functions.
 * \param meshfile Mesh file name.
 */
void CellData::readCellData(ifstream &meshfile) {
  readKeyword(meshfile);
  if (dims.get_ncell_data() > 0)
    readData(meshfile);
  readEndKeyword(meshfile);
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Reads and validates the cell_data block keyword.
 * \param meshfile Mesh file name.
 */
void CellData::readKeyword(ifstream &meshfile) {
  string dummyString;

  meshfile >> dummyString;
  Insist(dummyString == "celldat", "Invalid mesh file: celldat block missing");
  std::getline(meshfile, dummyString);
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Reads and validates the cell_data block data.
 * \param meshfile Mesh file name.
 */
void CellData::readData(ifstream &meshfile) {
  string dummyString;
  size_t cellNum;

  for (size_t i = 0; i < dims.get_ncells(); ++i) {
    meshfile >> cellNum;
    Insist(cellNum == i + 1, "Invalid mesh file: cell data index out of order");
    for (size_t j = 0; j < dims.get_ncell_data(); ++j) {
      meshfile >> data[i][j];
    }
    std::getline(meshfile, dummyString);
  }
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Reads and validate the end_celldat block keyword.
 * \param meshfile Mesh file name.
 */
void CellData::readEndKeyword(ifstream &meshfile) {
  string dummyString;

  meshfile >> dummyString;
  Insist(dummyString == "end_celldat", "Invalid mesh file: celldat block missing end");
  std::getline(meshfile, dummyString); // read and discard blank line.
}

} // end namespace rtt_RTT_Format_Reader

//------------------------------------------------------------------------------------------------//
// end of RTT_Format_Reader/CellData.cc
//------------------------------------------------------------------------------------------------//
