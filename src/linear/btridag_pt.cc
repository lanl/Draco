//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   linear/btridag_pt.cc
 * \author Kent Budge
 * \date   Wed Sep 15 13:03:41 MDT 2010
 * \brief  Implementation of block tridiagonal solver
 * \note   Copyright (C) 2016-2021 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "btridag.t.hh"
#include <vector>

namespace rtt_linear {
using namespace std;

template void btridag(vector<double> const &a, vector<double> const &b, vector<double> const &c,
                      vector<double> const &r, unsigned const n, unsigned const m,
                      vector<double> &u);

} // end namespace rtt_linear

//------------------------------------------------------------------------------------------------//
// end of btridag_pt.cc
//------------------------------------------------------------------------------------------------//
