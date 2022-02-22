//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   diagnostics/Diagnostics.cc
 * \author Thomas M. Evans
 * \date   Fri Dec  9 10:52:38 2005
 * \brief  Member definitions for Diagnostics class.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "Diagnostics.hh"
#include "ds++/Release.hh"
#include "ds++/fpe_trap.hh"
#include <sstream>

namespace rtt_diagnostics {

namespace Diagnostics {

//------------------------------------------------------------------------------------------------//
// DIAGNOSTICS NAMESPACE DEFINITIONS
//------------------------------------------------------------------------------------------------//

std::map<std::string, int> integers;
std::map<std::string, double> doubles;
std::map<std::string, std::vector<int>> vec_integers;
std::map<std::string, std::vector<double>> vec_doubles;

} // end namespace Diagnostics

} // end namespace rtt_diagnostics

//------------------------------------------------------------------------------------------------//
// end of Diagnostics.cc
//------------------------------------------------------------------------------------------------//
