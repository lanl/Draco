//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   linear/tred2_pt.cc
 * \author Kent Budge
 * \date   Thu Sep  2 14:49:55 2004
 * \brief  Specializations of tred2
 * \note   Copyright (C) 2016-2020 Triad National Security, LLC.
 *         All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "tred2.i.hh"
#include <cmath>
#include <vector>

namespace rtt_linear {

//------------------------------------------------------------------------------------------------//
// T1=T2=T3=vector<double>
//------------------------------------------------------------------------------------------------//
template void tred2(std::vector<double> &a, unsigned n, std::vector<double> &d,
                    std::vector<double> &e);

} // end namespace rtt_linear

//------------------------------------------------------------------------------------------------//
// end of tred2_pt.cc
//------------------------------------------------------------------------------------------------//
