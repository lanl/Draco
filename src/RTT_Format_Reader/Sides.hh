//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   RTT_Format_Reader/Sides.hh
 * \author B.T. Adams
 * \date   Wed Jun 7 10:33:26 2000
 * \brief  Header file for RTT_Format_Reader/Sides class.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_RTT_Format_Reader_Sides_hh
#define rtt_RTT_Format_Reader_Sides_hh

#include "CellDefs.hh"
#include "Nodes.hh"
#include "SideFlags.hh"

namespace rtt_RTT_Format_Reader {

//================================================================================================//
/*!
 * \class Sides
 * \brief Controls parsing, storing, and accessing the data specific to the sides block of the mesh
 *        file.
 */
//================================================================================================//
class Sides {

private:
  void readKeyword(std::ifstream &meshfile);
  void readData(std::ifstream &meshfile);
  void readEndKeyword(std::ifstream &meshfile);

  // typedefs
  using ifstream = std::ifstream;
  using string = std::string;
  using vector_int = std::vector<int>;
  using vector_vector_int = std::vector<std::vector<int>>;
  using vector_uint = std::vector<unsigned int>;
  using vector_vector_uint = std::vector<std::vector<unsigned int>>;

  const SideFlags &sideFlags;
  const Dims &dims;
  const CellDefs &cellDefs;
  vector_int sideType;
  vector_vector_uint nodes;
  vector_vector_int flags;

public:
  Sides(const SideFlags &sideFlags_, const Dims &dims_, const CellDefs &cellDefs_)
      : sideFlags(sideFlags_), dims(dims_), cellDefs(cellDefs_), sideType(dims.get_nsides()),
        nodes(dims.get_nsides()),
        flags(dims.get_nsides(), vector_int(dims.get_nside_flag_types())) { /* empty */
  }
  ~Sides() = default;
  Sides(Sides const &rhs) = delete;
  Sides(Sides &&rhs) noexcept = delete;
  Sides &operator=(Sides const &rhs) = delete;
  Sides &operator=(Sides &&rhs) noexcept = delete;

  void readSides(ifstream &meshfile);
  void redefineSides();

  /*!
   * \brief Returns the side type associated with the specified side.
   * \param side_numb Side number.
   * \return The side type.
   */
  int get_type(size_t side_numb) const { return sideType[side_numb]; }
  /*!
   * \brief Returns the node numbers associated with each side.
   * \return The node numbers for all of the sides.
   */
  vector_vector_uint get_nodes() const { return nodes; }
  /*!
   * \brief Returns the node numbers associated with the specified side.
   * \param side_numb Side number.
   * \return The side node numbers.
   */
  vector_uint get_nodes(size_t side_numb) const { return nodes[side_numb]; }
  /*!
   * \brief Returns the node number associated with the specified side and side-node index.
   * \param side_numb Side number.
   * \param node_numb Side-node index number.
   * \return The side node number.
   */
  int get_nodes(size_t side_numb, size_t node_numb) const { return nodes[side_numb][node_numb]; }
  /*!
   * \brief Returns the side flag for the specified side and flag index
   * \param side_numb Side number.
   * \param flag_numb Side flag index.
   * \return The side flag.
   */
  int get_flags(size_t side_numb, size_t flag_numb) const { return flags[side_numb][flag_numb]; }
};

} // end namespace rtt_RTT_Format_Reader

#endif // rtt_RTT_Format_Reader_Sides_hh

//------------------------------------------------------------------------------------------------//
// end of RTT_Format_Reader/Sides.hh
//------------------------------------------------------------------------------------------------//
