//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   RTT_Format_Reader/SideFlags.hh
 * \author B.T. Adams
 * \date   Wed Jun 7 10:33:26 2000
 * \brief  Header file for RTT_Format_Reader/SideFlags class.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_RTT_Format_Reader_SideFlags_hh
#define rtt_RTT_Format_Reader_SideFlags_hh

#include "Dims.hh"
#include "Flags.hh"
#include <memory>

namespace rtt_RTT_Format_Reader {

//================================================================================================//
/*!
 * \brief Controls parsing, storing, and accessing the data specific to the side flags block of the
 *        mesh file.
 */
//================================================================================================//
class SideFlags {

private:
  void readKeyword(std::ifstream &meshfile);
  void readFlagTypes(std::ifstream &meshfile);
  void readEndKeyword(std::ifstream &meshfile);

  // typedefs
  using ifstream = std::ifstream;
  using string = std::string;

  const Dims &dims;
  std::vector<std::shared_ptr<Flags>> flagTypes;

public:
  explicit SideFlags(const Dims &dims_) : dims(dims_), flagTypes(dims.get_nside_flag_types()) {}
  ~SideFlags() = default;
  SideFlags(SideFlags const &rhs) = delete;
  SideFlags(SideFlags &&rhs) noexcept = delete;
  SideFlags &operator=(SideFlags const &rhs) = delete;
  SideFlags &operator=(SideFlags &&rhs) noexcept = delete;

  void readSideFlags(ifstream &meshfile);

  /*!
   * \brief Validates the specified side flag type and number.
   * \param flagtype Side flag type number.
   * \param flag Flag number.
   * \return The existance of the side flag type and number.
   */
  bool allowed_flag(size_t flagtype, int flag) const {
    Insist(flagtype <= dims.get_nside_flag_types() - 1, "Invalid side flag type number!");
    return flagTypes[flagtype]->allowed_flag(flag);
  }

  /*!
   * \brief Returns the name of specified side flag type.
   * \param flagtype Side flag type number.
   * \return The side flag type name.
   */
  string get_flag_type(size_t flagtype) const {
    Insist(flagtype <= dims.get_nside_flag_types() - 1, "Invalid side flag type number!");
    return flagTypes[flagtype]->getFlagType();
  }

  int get_flag_type_index(string &desired_flag_type) const;

  /*!
   * \brief Returns the side flag number associated with the specified side flag type and side flag
   *        index.
   * \param flagtype Side flag type number.
   * \param flag_index Side flag index.
   * \return The side flag number.
   */
  int get_flag_number(size_t flagtype, size_t flag_index) const {
    Insist(flagtype <= dims.get_nside_flag_types() - 1, "Invalid side flag type number!");
    Insist(flag_index <= flagTypes[flagtype]->getFlagSize() - 1,
           "Invalid side flag number index number!");
    return flagTypes[flagtype]->getFlagNumber(flag_index);
  }

  /*!
   * \brief Returns the number of side flags for the specified side flag type.
   * \param flagtype Side flag type number.
   * \return The number of side flags.
   */
  size_t get_flag_size(size_t flagtype) const {
    Insist(flagtype <= dims.get_nside_flag_types() - 1, "Invalid side flag type number!");
    return flagTypes[flagtype]->getFlagSize();
  }
  /*!
   * \brief Returns the side flag name associated with the specified side flag index and side flag
   *        type.
   * \param flagtype Side flag type number.
   * \param flag_index Side flag index.
   * \return The side flag name.
   */
  string get_flag_name(size_t flagtype, size_t flag_index) const {
    Insist(flagtype <= dims.get_nside_flag_types() - 1, "Invalid side flag type number!");
    Insist(flag_index <= flagTypes[flagtype]->getFlagSize() - 1,
           "Invalid side flag name index number!");
    return flagTypes[flagtype]->getFlagName(flag_index);
  }
};

} // end namespace rtt_RTT_Format_Reader

#endif // rtt_RTT_Format_Reader_SideFlags_hh

//------------------------------------------------------------------------------------------------//
// end of RTT_Format_Reader/SideFlags.hh
//------------------------------------------------------------------------------------------------//
