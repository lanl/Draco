//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   RTT_Format_Reader/Nodes.hh
 * \author B.T. Adams
 * \date   Wed Jun 7 10:33:26 2000
 * \brief  Header file for RTT_Format_Reader/Nodes class.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_RTT_Format_Reader_Nodes_hh
#define rtt_RTT_Format_Reader_Nodes_hh

#include "NodeFlags.hh"

namespace rtt_RTT_Format_Reader {

//================================================================================================//
/*!
 * \brief Controls parsing, storing, and accessing the data specific to the nodes block of the mesh
 *        file.
 */
//================================================================================================//
class Nodes {

private:
  void readKeyword(std::ifstream &meshfile);
  void readData(std::ifstream &meshfile);
  void readEndKeyword(std::ifstream &meshfile);

  // typedefs
  using ifstream = std::ifstream;
  using string = std::string;
  using vector_int = std::vector<int>;
  using vector_vector_int = std::vector<std::vector<int>>;
  using vector_dbl = std::vector<double>;
  using vector_vector_dbl = std::vector<std::vector<double>>;

  const NodeFlags &nodeFlags;
  const Dims &dims;
  vector_vector_dbl coords;
  vector_int parents;
  vector_vector_int flags;

public:
  Nodes(const NodeFlags &nodeFlags_, const Dims &dims_)
      : nodeFlags(nodeFlags_), dims(dims_), coords(dims.get_nnodes(), vector_dbl(dims.get_ndim())),
        parents(dims.get_nnodes()),
        flags(dims.get_nnodes(), vector_int(dims.get_nnode_flag_types())) {}
  ~Nodes() = default;
  Nodes(Nodes const &rhs) = delete;
  Nodes(Nodes &&rhs) noexcept = delete;
  Nodes &operator=(Nodes const &rhs) = delete;
  Nodes &operator=(Nodes &&rhs) noexcept = delete;

  void readNodes(ifstream &meshfile);
  static int readNextInt(ifstream &meshfile);

  /*!
   * \brief Returns the coordinate values for each of the nodes.
   * \return The coordinate values for the nodes.
   */
  vector_vector_dbl get_coords() const { return coords; }

  /*!
   * \brief Returns all of the coordinate values for the specified node.
   * \param node_numb Node number.
   * \return The node coordinate values.
   */
  vector_dbl get_coords(size_t node_numb) const { return coords[node_numb]; }

  /*!
   * \brief Returns the coordinate value for the specified node and direction (i.e., x, y, and z).
   * \param node_numb Node number.
   * \param coord_index Coordinate index number (x = 0, y = 1, z = 2).
   * \return The node coordinate value.
   */
  double get_coords(size_t node_numb, size_t coord_index) const {
    return coords[node_numb][coord_index];
  }

  /*!
   * \brief Returns the node parent for the specified node.
   * \param node_numb Node number.
   * \return The node parent.
   */
  int get_parents(size_t node_numb) const { return parents[node_numb]; }

  /*!
   * \brief Returns the node flag for the specified node and flag index.
   * \param node_numb Node number.
   * \param flag_numb Node flag index.
   * \return The node flag.
   */
  int get_flags(size_t node_numb, size_t flag_numb) const { return flags[node_numb][flag_numb]; }
};

} // end namespace rtt_RTT_Format_Reader

#endif // rtt_RTT_Format_Reader_Nodes_hh

//------------------------------------------------------------------------------------------------//
// end of RTT_Format_Reader/Nodes.hh
//------------------------------------------------------------------------------------------------//
