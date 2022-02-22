//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   RTT_Format_Reader/NodeDataIDs.cc
 * \author B.T. Adams
 * \date   Wed Jun 7 10:33:26 2000
 * \brief  Implementation file for RTT_Format_Reader/NodeDataIDs class.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "NodeDataIDs.hh"

namespace rtt_RTT_Format_Reader {
/*!
 * \brief Parses the node_data_ids data block from the mesh file via calls to private member
 *        functions.
 * \param meshfile Mesh file name.
 */
void NodeDataIDs::readDataIDs(ifstream &meshfile) {
  readKeyword(meshfile);
  readData(meshfile);
  readEndKeyword(meshfile);
}
/*!
 * \brief Reads and validates the node_data_ids block keyword.
 * \param meshfile Mesh file name.
 */
void NodeDataIDs::readKeyword(ifstream &meshfile) {
  std::string dummyString;

  meshfile >> dummyString;
  Insist(dummyString == "node_data_ids", "Invalid mesh file: node_data_ids block missing");
  std::getline(meshfile, dummyString);
}
/*!
 * \brief Reads and validates the node_data_ids block data.
 * \param meshfile Mesh file name.
 */
void NodeDataIDs::readData(ifstream &meshfile) {
  int dataIDNum;
  string dummyString;

  for (size_t i = 0; i < dims.get_nnode_data(); ++i) {
    Check(i < names.size() && i < units.size());
    meshfile >> dataIDNum >> names[i] >> units[i];
    Insist(static_cast<size_t>(dataIDNum) == i + 1, "Invalid mesh file: node data ID out of order");
    std::getline(meshfile, dummyString);
  }
}
/*!
 * \brief Reads and validates the end_node_data_ids block keyword.
 * \param meshfile Mesh file name.
 */
void NodeDataIDs::readEndKeyword(ifstream &meshfile) {
  string dummyString;

  meshfile >> dummyString;
  Insist(dummyString == "end_node_data_ids", "Invalid mesh file: node_data_ids block missing end");
  std::getline(meshfile, dummyString); // read and discard blank line.
}

} // end namespace rtt_RTT_Format_Reader

//------------------------------------------------------------------------------------------------//
// end of RTT_Format_Reader/NodeDataIDs.cc
//------------------------------------------------------------------------------------------------//
