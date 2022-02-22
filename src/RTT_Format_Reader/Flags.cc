//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   RTT_Format_Reader/Flags.cc
 * \author B.T. Adams
 * \date   Wed Jun 7 10:33:26 2000
 * \brief  Implementation file for RTT_Format_Reader/Flags class.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "Flags.hh"

namespace rtt_RTT_Format_Reader {

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Used by the NodeFlags, SideFlags, and CellFlags class objects to parse the flag numbers
 *        and names.
 * \param meshfile Mesh file name.
 */
void Flags::readFlags(ifstream &meshfile) {
  string dummyString;

  for (size_t i = 0; i < nflags; ++i) {
    meshfile >> flag_nums[i] >> flag_names[i];
    std::getline(meshfile, dummyString);
  }
}

} // end namespace rtt_RTT_Format_Reader

//------------------------------------------------------------------------------------------------//
// end of RTT_Format_Reader/Flags.cc
//------------------------------------------------------------------------------------------------//
