//----------------------------------*-C++-*--------------------------------//
/*!
 * \file   RTT_Format_Reader/SideData.cc
 * \author B.T. Adams
 * \date   Wed Jun 7 10:33:26 2000
 * \brief  Implementation file for RTT_Format_Reader/SideData class.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "SideData.hh"

namespace rtt_RTT_Format_Reader {

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Parses the side_data block data from the mesh file via calls to private member functions.
 * \param meshfile Mesh file name.
 */
void SideData::readSideData(ifstream &meshfile) {
  readKeyword(meshfile);
  if (dims.get_nside_data() > 0)
    readData(meshfile);
  readEndKeyword(meshfile);
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Reads and validates the side_data block keyword.
 * \param meshfile Mesh file name.
 */
void SideData::readKeyword(ifstream &meshfile) {
  string dummyString;

  meshfile >> dummyString;
  Insist(dummyString == "sidedat", "Invalid mesh file: sidedat block missing");
  std::getline(meshfile, dummyString);
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Reads and validates the side data block data.
 * \param meshfile Mesh file name.
 */
void SideData::readData(ifstream &meshfile) {
  string dummyString;
  size_t sideNum;

  for (size_t i = 0; i < dims.get_nsides(); ++i) {
    meshfile >> sideNum;
    Insist(sideNum == i + 1, "Invalid mesh file: side data index out of order");
    for (size_t j = 0; j < dims.get_nside_data(); ++j)
      meshfile >> data[i][j];
    std::getline(meshfile, dummyString);
  }
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Reads and validates the end_sidedat block keyword.
 * \param meshfile Mesh file name.
 */
void SideData::readEndKeyword(ifstream &meshfile) {
  string dummyString;

  meshfile >> dummyString;
  Insist(dummyString == "end_sidedat", "Invalid mesh file: sidedat block missing end");
  std::getline(meshfile, dummyString); // read and discard blank line.
}

} // end namespace rtt_RTT_Format_Reader

//------------------------------------------------------------------------------------------------//
// end of RTT_Format_Reader/SideData.cc
//------------------------------------------------------------------------------------------------//
