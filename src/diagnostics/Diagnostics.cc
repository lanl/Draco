//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   diagnostics/Diagnostics.cc
 * \author Thomas M. Evans
 * \date   Fri Dec  9 10:52:38 2005
 * \brief  Member definitions for Diagnostics class.
 * \note   Copyright (C) 2010-2026 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "Diagnostics.hh"
#include <cstdint>
#include <map>
#include <string>
#include <vector>

namespace rtt_diagnostics::Diagnostics {

//------------------------------------------------------------------------------------------------//
// DIAGNOSTICS NAMESPACE DEFINITIONS
//------------------------------------------------------------------------------------------------//

std::map<std::string, int> integers;
std::map<std::string, int64_t> long_integers;
std::map<std::string, double> doubles;
std::map<std::string, std::vector<int>> vec_integers;
std::map<std::string, std::vector<double>> vec_doubles;

} // end namespace rtt_diagnostics::Diagnostics

//------------------------------------------------------------------------------------------------//
// end of Diagnostics.cc
//------------------------------------------------------------------------------------------------//
