//----------------------------------*-C++-*--------------------------------//
/*!
 * \file   RTT_Format_Reader/NodeData.cc
 * \author B.T. Adams
 * \date   Wed Jun 7 10:33:26 2000
 * \brief  Implementation file for RTT_Format_Reader/NodeData class.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "NodeData.hh"

namespace rtt_RTT_Format_Reader {

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Parses the node_data block data from the mesh file via calls to private member functions.
 * \param meshfile Mesh file name.
 */
void NodeData::readNodeData(ifstream &meshfile) {
  readKeyword(meshfile);
  if (dims.get_nnode_data() > 0)
    readData(meshfile);
  readEndKeyword(meshfile);
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Reads and validates the node_data block keyword.
 * \param meshfile Mesh file name.
 */
void NodeData::readKeyword(ifstream &meshfile) {
  string dummyString;

  meshfile >> dummyString;
  Insist(dummyString == "nodedat", "Invalid mesh file: nodedat block missing");
  std::getline(meshfile, dummyString);
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Reads and validates the node_data block datae.
 * \param meshfile Mesh file name.
 */
void NodeData::readData(ifstream &meshfile) {
  string dummyString;
  size_t nodeNum;

  for (size_t i = 0; i < dims.get_nnodes(); ++i) {
    meshfile >> nodeNum;
    Insist(nodeNum == i + 1, "Invalid mesh file: node data index out of order");
    for (size_t j = 0; j < dims.get_nnode_data(); ++j)
      meshfile >> data[i][j];
    std::getline(meshfile, dummyString);
  }
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Reads and validates the end_nodedat block keyword.
 * \param meshfile Mesh file name.
 */
void NodeData::readEndKeyword(ifstream &meshfile) {
  string dummyString;

  meshfile >> dummyString;
  Insist(dummyString == "end_nodedat", "Invalid mesh file: nodedat block missing end");
  std::getline(meshfile, dummyString); // read and discard blank line.
}

} // end namespace rtt_RTT_Format_Reader

//------------------------------------------------------------------------------------------------//
// end of RTT_Format_Reader/NodeData.cc
//------------------------------------------------------------------------------------------------//
