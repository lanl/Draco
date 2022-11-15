//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   c4/C4_MPI_get_processor_name.cc
 * \author Thomas M. Evans
 * \date   Thu Mar 21 16:56:17 2002
 * \brief  C4 MPI implementation.
 * \note   Copyright (C) 2014-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "c4/config.h"
#include <array>
#include <string>

#ifdef C4_MPI

#include "C4_Functions.hh"

namespace rtt_c4 {
//------------------------------------------------------------------------------------------------//
// get_processor_name
//------------------------------------------------------------------------------------------------//
std::string get_processor_name() {
  int namelen(0);
  // char processor_name[DRACO_MAX_PROCESSOR_NAME];
  std::array<char, DRACO_MAX_PROCESSOR_NAME> processor_name{};
  MPI_Get_processor_name(processor_name.data(), &namelen);
  std::string pname(processor_name.begin(), processor_name.begin() + namelen);
  Ensure(pname.size() == static_cast<size_t>(namelen));
  return pname;
}

} // end namespace rtt_c4

#endif // C4_MPI

//------------------------------------------------------------------------------------------------//
// end of C4_MPI_get_processor_name.cc
//------------------------------------------------------------------------------------------------//
