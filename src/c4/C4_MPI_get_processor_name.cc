//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   c4/C4_MPI_get_processor_name.cc
 * \author Thomas M. Evans
 * \brief  C4 MPI implementation.
 * \note   Copyright (C) 2014-2026 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "c4/config.h"

#ifdef C4_MPI

#include "dsxx/Assert.hh"
#include "c4/C4_Functions.hh"
#include <array>
#include <cstddef> // size_t
#include <string>

namespace rtt_c4 {
//------------------------------------------------------------------------------------------------//
// get_processor_name
//------------------------------------------------------------------------------------------------//
std::string get_processor_name() {
  int namelen(0);
  // char processor_name[DRACO_MAX_PROCESSOR_NAME];
  std::array<char, DRACO_MAX_PROCESSOR_NAME> processor_name{}; // NOLINT(misc-include-cleaner)
  MPI_Get_processor_name(processor_name.data(), &namelen);     // NOLINT(misc-include-cleaner)
  std::string pname(processor_name.begin(), processor_name.begin() + namelen);
  Ensure(pname.size() == static_cast<size_t>(namelen));
  return pname;
}

} // end namespace rtt_c4

#endif // C4_MPI

//------------------------------------------------------------------------------------------------//
// end of C4_MPI_get_processor_name.cc
//------------------------------------------------------------------------------------------------//
