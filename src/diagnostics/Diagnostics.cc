//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   diagnostics/Diagnostics.cc
 * \author Thomas M. Evans
 * \date   Fri Dec  9 10:52:38 2005
 * \brief  Member definitions for Diagnostics class.
 * \note   Copyright (C) 2010-2024 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "Diagnostics.hh"
#include "dsxx/Release.hh"
#include "dsxx/fpe_trap.hh"
#include <sstream>

#if defined(__cplusplus) && __cplusplus > 201402L // C++17
namespace rtt_diagnostics::Diagnostics {
#else
namespace rtt_diagnostics {
namespace Diagnostics {
#endif

//------------------------------------------------------------------------------------------------//
// DIAGNOSTICS NAMESPACE DEFINITIONS
//------------------------------------------------------------------------------------------------//

std::map<std::string, int> integers;
std::map<std::string, double> doubles;
std::map<std::string, std::vector<int>> vec_integers;
std::map<std::string, std::vector<double>> vec_doubles;

#if defined(__cplusplus) && __cplusplus > 201402L // C++17
} // end namespace rtt_diagnostics::Diagnostics
#else
} // end namespace Diagnostics
} // end namespace rtt_diagnostics
#endif

//------------------------------------------------------------------------------------------------//
// end of Diagnostics.cc
//------------------------------------------------------------------------------------------------//
