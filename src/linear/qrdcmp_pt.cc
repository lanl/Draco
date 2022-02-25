//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   linear/qrdcmp_pt.cc
 * \author Kent Budge
 * \date   Wed Aug 11 15:21:38 2004
 * \brief  Specializations of qrdcmp
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "qrdcmp.t.hh"
#include <vector>

namespace rtt_linear {

//------------------------------------------------------------------------------------------------//
// RandomContainer=vector<double>
//------------------------------------------------------------------------------------------------//

template bool qrdcmp(std::vector<double> &a, unsigned n, std::vector<double> &c,
                     std::vector<double> &d);

} // end namespace rtt_linear

//------------------------------------------------------------------------------------------------//
// end of qrdcmp_pt.cc
//------------------------------------------------------------------------------------------------//
