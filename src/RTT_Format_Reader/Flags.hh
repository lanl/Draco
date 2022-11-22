//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   RTT_Format_Reader/Flags.hh
 * \author B.T. Adams
 * \date   Wed Jun 7 10:33:26 2000
 * \brief  Header file for RTT_Format_Reader/Flags class.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_RTT_Format_Reader_Flags_hh
#define rtt_RTT_Format_Reader_Flags_hh

#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

namespace rtt_RTT_Format_Reader {

//================================================================================================//
/*!
 * \class Flags
 * \brief Controls parsing, storing, and accessing the data contained in the node, side, and cell
 *        flag blocks of the mesh file.
 */
//================================================================================================//
class Flags {

  // typedefs
  using ifstream = std::ifstream;
  using string = std::string;
  using vector_str = std::vector<string>;
  using vector_int = std::vector<int>;

  size_t nflags;
  string name;
  vector_int flag_nums;
  vector_str flag_names;

public:
  Flags(size_t nflags_, string name_)
      : nflags(nflags_), name(std::move(name_)), flag_nums(nflags), flag_names(nflags) {}
  ~Flags() = default;
  Flags(Flags const &rhs) = delete;
  Flags(Flags &&rhs) noexcept = delete;
  Flags &operator=(Flags const &rhs) = delete;
  Flags &operator=(Flags &&rhs) noexcept = delete;

  void readFlags(ifstream &meshfile);

  /*!
   * \brief Validates the specified flag index.
   * \param flag Flag index.
   * \return The existence of the flag.
   */
  bool allowed_flag(size_t flag) const {
    return flag_nums.end() != std::find(flag_nums.begin(), flag_nums.end(), static_cast<int>(flag));
  }

  /*!
   * \brief Returns the name of specified flag type
   * \return The Flag type name.
   */
  string getFlagType() const { return name; }

  /*!
   * \brief Returns the flag number associated with the specified flag index.
   * \param flag Flag index.
   * \return The Flag number.
   */
  int getFlagNumber(size_t flag) const { return flag_nums[flag]; }

  /*!
   * \brief Returns the flag name associated with the specified index.
   * \param flag Flag index.
   * \return The flag name.
   */
  string getFlagName(size_t flag) const { return flag_names[flag]; }

  /*!
   * \brief Returns the number of flags.
   * \return The number of flags.
   */
  size_t getFlagSize() const { return nflags; }
};

} // end namespace rtt_RTT_Format_Reader

#endif // rtt_RTT_Format_Reader_Flags_hh

//------------------------------------------------------------------------------------------------//
// end of RTT_Format_Reader/Flags.hh
//------------------------------------------------------------------------------------------------//
