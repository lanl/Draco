//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   linear/qrupdt_pt.cc
 * \author Kent Budge
 * \date   Wed Aug 11 15:21:38 2004
 * \brief  Specializations of qrupdt
 * \note   Copyright (C) 2004-2021 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "qrupdt.t.hh"
#include <vector>

namespace rtt_linear {
using std::vector;

//------------------------------------------------------------------------------------------------//
// T = RandomContainer = vector<double>
//------------------------------------------------------------------------------------------------//

template void qrupdt(vector<double> &r, vector<double> &qt, const unsigned n, vector<double> &u,
                     vector<double> &v);

} // end namespace rtt_linear

//------------------------------------------------------------------------------------------------//
//  end of linear/qrupdt_pt.cc
//------------------------------------------------------------------------------------------------//
