//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   diagnostics/Diagnostics.cc
 * \author Thomas M. Evans
 * \date   Fri Dec  9 10:52:38 2005
 * \brief  Member definitions for Diagnostics class.
 * \note   Copyright (C) 2010-2024 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "Diagnostics.hh"
#include "ds++/Release.hh"
#include "ds++/fpe_trap.hh"
#include <sstream>

#ifdef __cplusplus
#if __cplusplus > 201402L
namespace rtt_diagnostics::Diagnostics {
#else
namespace rtt_diagnostics {
namespace Diagnostics {
#endif
#endif

//------------------------------------------------------------------------------------------------//
// DIAGNOSTICS NAMESPACE DEFINITIONS
//------------------------------------------------------------------------------------------------//

std::map<std::string, int> integers;
std::map<std::string, double> doubles;
std::map<std::string, std::vector<int>> vec_integers;
std::map<std::string, std::vector<double>> vec_doubles;

#ifdef __cplusplus
#if __cplusplus > 201402L
} // end namespace rtt_diagnostics::Diagnostics
#else
} // end namespace Diagnostics
} // end namespace rtt_diagnostics
#endif
#endif

//------------------------------------------------------------------------------------------------//
// end of Diagnostics.cc
//------------------------------------------------------------------------------------------------//
