//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   RTT_Format_Reader/SideFlags.hh
 * \author B.T. Adams
 * \date   Wed Jun 7 10:33:26 2000
 * \brief  Header file for RTT_Format_Reader/SideFlags class.
 * \note   Copyright (C) 2016-2017 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __RTT_Format_Reader_SideFlags_hh__
#define __RTT_Format_Reader_SideFlags_hh__

#include "Dims.hh"
#include "Flags.hh"
#include "ds++/SP.hh"

namespace rtt_RTT_Format_Reader {
/*!
 * \brief Controls parsing, storing, and accessing the data specific to the
 *        side flags block of the mesh file.
 */
class SideFlags {
  // typedefs
  typedef std::ifstream ifstream;
  typedef std::string string;

  const Dims &dims;
  std::vector<rtt_dsxx::SP<Flags>> flagTypes;

public:
  SideFlags(const Dims &dims_)
      : dims(dims_), flagTypes(dims.get_nside_flag_types()) {}
  ~SideFlags() {}

  void readSideFlags(ifstream &meshfile);

private:
  void readKeyword(ifstream &meshfile);
  void readFlagTypes(ifstream &meshfile);
  void readEndKeyword(ifstream &meshfile);

public:
  /*!
 * \brief Validates the specified side flag type and number.
 * \param flagtype Side flag type number.
 * \param flag Flag number.
 * \return The existance of the side flag type and number.
 */
  bool allowed_flag(int flagtype, int flag) const {
    Insist(flagtype <= dims.get_nside_flag_types() - 1,
           "Invalid side flag type number!");
    return flagTypes[flagtype]->allowed_flag(flag);
  }
  /*!
 * \brief Returns the name of specified side flag type.
 * \param flagtype Side flag type number.
 * \return The side flag type name.
 */
  string get_flag_type(int flagtype) const {
    Insist(flagtype <= dims.get_nside_flag_types() - 1,
           "Invalid side flag type number!");
    return flagTypes[flagtype]->getFlagType();
  }

  int get_flag_type_index(string &desired_flag_type) const;
  /*!
 * \brief Returns the side flag number associated with the specified side flag
 *        type and side flag index.
 * \param flagtype Side flag type number.
 * \param flag_index Side flag index.
 * \return The side flag number.
 */
  int get_flag_number(int flagtype, int flag_index) const {
    Insist(flagtype <= dims.get_nside_flag_types() - 1,
           "Invalid side flag type number!");
    Insist(flag_index <= flagTypes[flagtype]->getFlagSize() - 1,
           "Invalid side flag number index number!");
    return flagTypes[flagtype]->getFlagNumber(flag_index);
  }
  /*!
 * \brief Returns the number of side flags for the specified side flag type.
 * \param flagtype Side flag type number.
 * \return The number of side flags.
 */
  int get_flag_size(int flagtype) const {
    Insist(flagtype <= dims.get_nside_flag_types() - 1,
           "Invalid side flag type number!");
    return flagTypes[flagtype]->getFlagSize();
  }
  /*!
 * \brief Returns the side flag name associated with the specified side flag
 *        index and side flag type.
 * \param flagtype Side flag type number.
 * \param flag_index Side flag index.
 * \return The side flag name.
 */
  string get_flag_name(int flagtype, int flag_index) const {
    Insist(flagtype <= dims.get_nside_flag_types() - 1,
           "Invalid side flag type number!");
    Insist(flag_index <= flagTypes[flagtype]->getFlagSize() - 1,
           "Invalid side flag name index number!");
    return flagTypes[flagtype]->getFlagName(flag_index);
  }
};

} // end namespace rtt_RTT_Format_Reader

#endif // __RTT_Format_Reader_SideFlags_hh__

//---------------------------------------------------------------------------//
// end of RTT_Format_Reader/SideFlags.hh
//---------------------------------------------------------------------------//
