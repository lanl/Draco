//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   linear/gaussj_pt.cc
 * \author Kent Budge
 * \date   Wed Aug 11 15:21:38 2004
 * \brief  Specializations of gaussj
 * \note   Copyright (C) 2016-2021 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "gaussj.t.hh"
#include <vector>

namespace rtt_linear {
using std::vector;

template void gaussj(vector<double> &A, unsigned n, vector<double> &b, unsigned m);
template void gaussj(vector<vector<double>> &A, vector<double> &b);

} // end namespace rtt_linear

//------------------------------------------------------------------------------------------------//
// end of linear/gaussj_pt.cc
//------------------------------------------------------------------------------------------------//
