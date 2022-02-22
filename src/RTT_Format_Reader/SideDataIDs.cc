//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   RTT_Format_Reader/SideDataIDs.cc
 * \author B.T. Adams
 * \date   Wed Jun 7 10:33:26 2000
 * \brief  Implementation file for RTT_Format_Reader/SideDataIDs class.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "SideDataIDs.hh"

namespace rtt_RTT_Format_Reader {
/*!
 * \brief Parses the side_data_ids data block from the mesh file via calls to private member
 *        functions.
 * \param meshfile Mesh file name.
 */
void SideDataIDs::readDataIDs(ifstream &meshfile) {
  readKeyword(meshfile);
  readData(meshfile);
  readEndKeyword(meshfile);
}
/*!
 * \brief Reads and validates the side_data_ids block keyword.
 * \param meshfile Mesh file name.
 */
void SideDataIDs::readKeyword(ifstream &meshfile) {
  string dummyString;

  meshfile >> dummyString;
  Insist(dummyString == "side_data_ids", "Invalid mesh file: side_data_ids block missing");
  std::getline(meshfile, dummyString);
}
/*!
 * \brief Reads and validates the side_data_ids block data.
 * \param meshfile Mesh file name.
 */
void SideDataIDs::readData(ifstream &meshfile) {
  int dataIDNum;
  string dummyString;

  for (unsigned i = 0; i < static_cast<unsigned int>(dims.get_nside_data()); ++i) {
    Check(i < names.size() && i < units.size());
    meshfile >> dataIDNum >> names[i] >> units[i];
    Insist(static_cast<unsigned int>(dataIDNum) == i + 1,
           "Invalid mesh file: side data ID out of order");
    std::getline(meshfile, dummyString);
  }
}
/*!
 * \brief Reads and validates the end_side_data_ids block keyword.
 * \param meshfile Mesh file name.
 */
void SideDataIDs::readEndKeyword(ifstream &meshfile) {
  string dummyString;

  meshfile >> dummyString;
  Insist(dummyString == "end_side_data_ids", "Invalid mesh file: side_data_ids block missing end");
  std::getline(meshfile, dummyString); // read and discard blank line.
}

} // end namespace rtt_RTT_Format_Reader

//------------------------------------------------------------------------------------------------//
// end of RTT_Format_Reader/SideDataIDs.cc
//------------------------------------------------------------------------------------------------//
