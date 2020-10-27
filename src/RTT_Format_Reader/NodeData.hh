//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   RTT_Format_Reader/NodeData.hh
 * \author B.T. Adams
 * \date   Wed Jun 7 10:33:26 2000
 * \brief  Header file for RTT_Format_Reader/NodeData class.
 * \note   Copyright (C) 2016-2020 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_RTT_Format_Reader_NodeData_hh
#define rtt_RTT_Format_Reader_NodeData_hh

#include "Nodes.hh"

namespace rtt_RTT_Format_Reader {
/*!
 * \brief Controls parsing, storing, and accessing the data specific to the nodedata block of the
 *        mesh file.
 */
class NodeData {
  // typedefs
  typedef std::ifstream ifstream;
  typedef std::string string;
  typedef std::vector<double> vector_dbl;
  typedef std::vector<std::vector<double>> vector_vector_dbl;

  const Dims &dims;
  vector_vector_dbl data;

public:
  NodeData(const Dims &dims_)
      : dims(dims_), data(dims.get_nnodes(), vector_dbl(dims.get_nnode_data())) {}
  ~NodeData() {}

  void readNodeData(ifstream &meshfile);

private:
  void readKeyword(ifstream &meshfile);
  void readData(ifstream &meshfile);
  void readEndKeyword(ifstream &meshfile);

public:
  /*!
   * \brief Returns all of the data field values for each of the nodes.
   * \return The data field values for each of the nodes.
   */
  vector_vector_dbl get_data() const { return data; }

  /*!
   * \brief Returns all of the data field values for the specified node.
   * \param node_numb Node number.
   * \return The node data field values.
   */
  vector_dbl get_data(size_t node_numb) const { return data[node_numb]; }

  /*!
   * \brief Returns the specified data field value for the specified node.
   * \param node_numb Node number.
   * \param data_index Data field.
   * \return The node data field value.
   */
  double get_data(size_t node_numb, size_t data_index) const { return data[node_numb][data_index]; }
};

} // end namespace rtt_RTT_Format_Reader

#endif // rtt_RTT_Format_Reader_NodeData_hh

//------------------------------------------------------------------------------------------------//
// end of RTT_Format_Reader/NodeData.hh
//------------------------------------------------------------------------------------------------//
