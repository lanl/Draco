//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   RTT_Format_Reader/SideDataIDs.hh
 * \author Shawn Pautz/B.T. Adams
 * \date   Wed Jun 7 10:33:26 2000
 * \brief  Header file for RTT_Format_Reader/SideDataIDs class.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_RTT_Format_Reader_SideDataIDs_hh
#define rtt_RTT_Format_Reader_SideDataIDs_hh

#include "Dims.hh"
#include "ds++/Assert.hh"
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

namespace rtt_RTT_Format_Reader {

//================================================================================================//
/*!
 * \brief Controls parsing, storing, and accessing the data specific to the side data ids block of
 *        the mesh file.
 */
//================================================================================================//
class SideDataIDs {

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
  explicit SideDataIDs(const Dims &dims_)
      : dims(dims_), names(dims.get_nside_data()), units(dims.get_nside_data()) {}
  ~SideDataIDs() = default;
  SideDataIDs(SideDataIDs const &rhs) = delete;
  SideDataIDs(SideDataIDs &&rhs) noexcept = delete;
  SideDataIDs &operator=(SideDataIDs const &rhs) = delete;
  SideDataIDs &operator=(SideDataIDs &&rhs) noexcept = delete;

  void readDataIDs(ifstream &meshfile);

  /*!
   * \brief Returns the specified side_data_id name.
   * \param id_numb side_data_id index number.
   * \return The side_data_id name.
   */
  string get_data_id_name(size_t id_numb) const {
    Insist(id_numb <= dims.get_nside_data() - 1, "Invalid side data id number!");
    return names[id_numb];
  }

  /*!
   * \brief Returns the units associated with the specified side_data_id.
   * \param id_numb side_data_id index number.
   * \return The side_data_id units.
   */
  string get_data_id_units(size_t id_numb) const {
    Insist(id_numb <= dims.get_nside_data() - 1, "Invalid side data id number!");
    return units[id_numb];
  }
};

} // end namespace rtt_RTT_Format_Reader

#endif // rtt_RTT_Format_Reader_SideDataIDs_hh

//------------------------------------------------------------------------------------------------//
// end of RTT_Format_Reader/SideDataIDs.hh
//------------------------------------------------------------------------------------------------//
