//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   RTT_Format_Reader/NodeFlags.cc
 * \author B.T. Adams
 * \date   Wed Jun 7 10:33:26 2000
 * \brief  Implementation file for RTT_Format_Reader/NodeFlags class.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "NodeFlags.hh"

namespace rtt_RTT_Format_Reader {

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Parses the node_flags data block from the mesh file via calls to private member functions.
 * \param meshfile Mesh file name.
 */
void NodeFlags::readNodeFlags(ifstream &meshfile) {
  readKeyword(meshfile);
  readFlagTypes(meshfile);
  readEndKeyword(meshfile);
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Reads and validates the node_flags block keyword.
 * \param meshfile Mesh file name.
 */
void NodeFlags::readKeyword(ifstream &meshfile) {
  string dummyString;

  meshfile >> dummyString;
  Insist(dummyString == "node_flags", "Invalid mesh file: node_flags block missing");
  std::getline(meshfile, dummyString);
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Reads and validates the node_flags block data.
 * \param meshfile Mesh file name.
 */
void NodeFlags::readFlagTypes(ifstream &meshfile) {
  int flagTypeNum;
  string dummyString;

  for (size_t i = 0; i < dims.get_nnode_flag_types(); ++i) {
    meshfile >> flagTypeNum >> dummyString;
    Insist(static_cast<size_t>(flagTypeNum) == i + 1,
           "Invalid mesh file: node flag type out of order");
    Check(i < flagTypes.size());
    flagTypes[i] = std::make_shared<Flags>(dims.get_nnode_flags(i), dummyString);
    std::getline(meshfile, dummyString);
    flagTypes[i]->readFlags(meshfile);
  }
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Reads and validates the end_node_flags block keyword.
 * \param meshfile Mesh file name.
 */
void NodeFlags::readEndKeyword(ifstream &meshfile) {
  string dummyString;

  meshfile >> dummyString;
  Insist(dummyString == "end_node_flags", "Invalid mesh file: node_flags block missing end");
  std::getline(meshfile, dummyString); // read and discard blank line.
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Returns the index to the node flag type that contains the specified string.
 * \param desired_flag_type Flag type.
 * \return The node flag type index.
 */
int NodeFlags::get_flag_type_index(string &desired_flag_type) const {
  int flag_type_index = -1;
  for (size_t f = 0; f < dims.get_nnode_flag_types(); f++) {
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
// end of RTT_Format_Reader/NodeFlags.cc
//------------------------------------------------------------------------------------------------//
