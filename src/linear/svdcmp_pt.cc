//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   linear/svdcmp_pt.cc
 * \author Kent Budge
 * \date   Wed Aug 11 15:21:38 2004
 * \brief  Specializations of svdcmp
 * \note   Copyright (C) 2016-2021 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "svdcmp.t.hh"
#include <vector>

namespace rtt_linear {
using std::vector;

template void svdcmp(vector<double> &a, const unsigned m, const unsigned n, vector<double> &w,
                     vector<double> &v);

} // end namespace rtt_linear

//------------------------------------------------------------------------------------------------//
// end of svdcmp.cc
//------------------------------------------------------------------------------------------------//
