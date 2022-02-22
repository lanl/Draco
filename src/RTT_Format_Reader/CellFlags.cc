//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   RTT_Format_Reader/CellFlags.cc
 * \author B.T. Adams
 * \date   Mon Jun 7 10:33:26 2000
 * \brief  Implementation file for RTT_Format_Reader/CellFlags class
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "CellFlags.hh"

namespace rtt_RTT_Format_Reader {

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Parses the cell_flags data block of the mesh file via calls to private member functions.
 * \param meshfile Mesh file name.
 */
void CellFlags::readCellFlags(ifstream &meshfile) {
  readKeyword(meshfile);
  readFlagTypes(meshfile);
  readEndKeyword(meshfile);
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Reads and validates the cell_flags block keyword.
 * \param meshfile Mesh file name.
 */
void CellFlags::readKeyword(ifstream &meshfile) {
  string dummyString;

  meshfile >> dummyString;
  Insist(dummyString == "cell_flags", "Invalid mesh file: cell_flags block missing");
  std::getline(meshfile, dummyString);
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Reads and validates the cell_flags block data.
 * \param meshfile Mesh file name.
 */
void CellFlags::readFlagTypes(ifstream &meshfile) {
  int flagTypeNum;
  string dummyString;

  for (size_t i = 0; i < dims.get_ncell_flag_types(); ++i) {
    meshfile >> flagTypeNum >> dummyString;
    Insist(static_cast<size_t>(flagTypeNum) == i + 1,
           "Invalid mesh file: cell flag type out of order");
    Check(i < flagTypes.size());
    Check(i < INT_MAX);
    flagTypes[i] = std::make_shared<Flags>(dims.get_ncell_flags(static_cast<int>(i)), dummyString);
    std::getline(meshfile, dummyString);
    flagTypes[i]->readFlags(meshfile);
  }
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Reads and validates the end_cell_flags block keyword.
 * \param meshfile Mesh file name.
 */
void CellFlags::readEndKeyword(ifstream &meshfile) {
  string dummyString;

  meshfile >> dummyString;
  Insist(dummyString == "end_cell_flags", "Invalid mesh file: cell_flags block missing end");
  std::getline(meshfile, dummyString); // read and discard blank line.
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Returns the index to the cell flag type that contains the specified string.
 * \param desired_flag_type Flag type.
 * \return The cell flag type index.
 */
int CellFlags::get_flag_type_index(string &desired_flag_type) const {
  int flag_type_index = -1;
  for (size_t f = 0; f < dims.get_ncell_flag_types(); f++) {
    string flag_type = flagTypes[f]->getFlagType();
    if (flag_type == desired_flag_type) {
      Check(f < INT_MAX);
      flag_type_index = static_cast<int>(f);
    }
  }
  return flag_type_index;
}

} // end namespace rtt_RTT_Format_Reader

//------------------------------------------------------------------------------------------------//
// end of RTT_Format_Reader/CellFlags.cc
//------------------------------------------------------------------------------------------------//
