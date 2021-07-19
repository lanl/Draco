//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   linear/svbksb_pt.cc
 * \author Kent Budge
 * \date   Tue Aug 10 13:08:03 2004
 * \brief  Specializations of svbksb
 * \note   Copyright (C) 2016-2021 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "svbksb.t.hh"
#include <vector>

namespace rtt_linear {
using std::vector;

//------------------------------------------------------------------------------------------------//
// T=vector<double>
//------------------------------------------------------------------------------------------------//

template void svbksb(const vector<double> &u, const vector<double> &w, const vector<double> &v,
                     const unsigned m, const unsigned n, const vector<double> &b,
                     vector<double> &x);

} // end namespace rtt_linear

//------------------------------------------------------------------------------------------------//
// end of svbksb_pt.cc
//------------------------------------------------------------------------------------------------//
