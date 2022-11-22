//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   RTT_Format_Reader/NodeDataIDs.hh
 * \author B.T. Adams
 * \date   Wed Jun 7 10:33:26 2000
 * \brief  Header file for RTT_Format_Reader/NodeDataIDs class.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_RTT_Format_Reader_NodeDataIDs_hh
#define rtt_RTT_Format_Reader_NodeDataIDs_hh

#include "Dims.hh"

namespace rtt_RTT_Format_Reader {

//================================================================================================//
/*!
 * \brief Controls parsing, storing, and accessing the data specific to the node data ids block of
 *        the mesh file.
 */
//================================================================================================//
class NodeDataIDs {

private:
  void readKeyword(std::ifstream &meshfile);
  void readData(std::ifstream &meshfile);
  void readEndKeyword(std::ifstream &meshfile);

  // typedefs
  using ifstream = std::ifstream;
  using string = std::string;
  using vector_str = std::vector<string>;

  const Dims &dims;
  vector_str names;
  vector_str units;

public:
  explicit NodeDataIDs(const Dims &dims_)
      : dims(dims_), names(dims.get_nnode_data()), units(dims.get_nnode_data()) {}
  ~NodeDataIDs() = default;
  NodeDataIDs(NodeDataIDs const &rhs) = delete;
  NodeDataIDs(NodeDataIDs &&rhs) noexcept = delete;
  NodeDataIDs &operator=(NodeDataIDs const &rhs) = delete;
  NodeDataIDs &operator=(NodeDataIDs &&rhs) noexcept = delete;

  void readDataIDs(ifstream &meshfile);

  /*!
   * \brief Returns the specified node_data_id name.
   * \param id_numb node_data_id index number.
   * \return The node_data_id name.
   */
  string get_data_id_name(size_t id_numb) const {
    Insist(id_numb <= dims.get_nnode_data() - 1, "Invalid node data id number!");
    return names[id_numb];
  }

  /*!
   * \brief Returns the units associated with the specified node_data_id.
   * \param id_numb node_data_id index number.
   * \return The node_data_id units.
   */
  string get_data_id_units(size_t id_numb) const {
    Insist(id_numb <= dims.get_nnode_data() - 1, "Invalid node data id number!");
    return units[id_numb];
  }
};

} // end namespace rtt_RTT_Format_Reader

#endif // rtt_RTT_Format_Reader_NodeDataIDs_hh

//------------------------------------------------------------------------------------------------//
// end of RTT_Format_Reader/NodeDataIDs.hh
//------------------------------------------------------------------------------------------------//
