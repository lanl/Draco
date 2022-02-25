//--------------------------------------------*-C++-*---------------------------------------------//
/*! \file   UnitSystemEnums.cc
 *  \author Kelly Thompson
 *  \brief  This file contains enums, conversion factors and labels that help define a UnitSystem.
 *  \date   Mon Nov 03 20:54:05 2003
 *  \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "UnitSystemEnums.hh"
#include "ds++/Assert.hh"
#include <string>
#include <vector>

namespace rtt_units {

//------------------------------------------------------------------------------------------------//
// HELPER FUNCTIONS
//------------------------------------------------------------------------------------------------//

/*!
 * \brief Extract unit labels from list in UnitSystemEnums.hh.
 *
 * The unit labels are stored in UnitSystemEnum.hh as a std::string.  Each label is dilemited with a
 * comma.
 *
 * \param[in] pos Provide the enum value for the current unit.  pos is the integer equivalent of the
 *              enum value.
 * \param[in] labels The string that contains all of the label names.  Each name dilemited with a
 *              comma.
 * \return A std::string that contains only one label
 */
std::string setUnitLabel(size_t const pos, std::string const &labels) {

  Require(labels.length() > 0);

  // Store the location of the first letter of each label.  Also append the position for one past
  // the end of the original string.
  std::vector<std::string::size_type> word_positions;

  // idx is the index for the first character of the lable.  numChars is the length of the label.
  std::string::size_type idx(0), numChars(0);
  word_positions.push_back(0);
  while ((idx = labels.find(',', idx)) != std::string::npos)
    word_positions.push_back(++idx);

  // append one past the end the string.  This lets us compute the label length without using if
  // tests.
  word_positions.push_back(labels.length() + 1);
  idx = word_positions[pos];
  numChars = word_positions[pos + 1] - idx - 1;
  std::string retVal = labels.substr(idx, numChars);

  return retVal;

} // end setUnitLabel

} // end namespace rtt_units

//------------------------------------------------------------------------------------------------//
// end of UnitSystemEnums.cc
//------------------------------------------------------------------------------------------------//
