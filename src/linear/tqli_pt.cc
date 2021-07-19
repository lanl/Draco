//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   linear/tqli_pt.cc
 * \author Kent Budge
 * \date   Thu Sep  2 15:00:32 2004
 * \brief  Specializations of tqli
 * \note   Copyright (C) 2016-2021 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "tqli.t.hh"
#include <vector>

namespace rtt_linear {
using std::vector;

//------------------------------------------------------------------------------------------------//
// T=vector<double>
//------------------------------------------------------------------------------------------------//

template void tqli(vector<double> &d, vector<double> &e, const unsigned n, vector<double> &z);

} // end namespace rtt_linear

//------------------------------------------------------------------------------------------------//
// end of tqli.cc
//------------------------------------------------------------------------------------------------//
