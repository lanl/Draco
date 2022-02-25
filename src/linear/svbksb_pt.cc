//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   linear/svbksb_pt.cc
 * \author Kent Budge
 * \date   Tue Aug 10 13:08:03 2004
 * \brief  Specializations of svbksb
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "svbksb.t.hh"
#include <vector>

namespace rtt_linear {

//------------------------------------------------------------------------------------------------//
// T=vector<double>
//------------------------------------------------------------------------------------------------//

template void svbksb(const std::vector<double> &u, const std::vector<double> &w,
                     const std::vector<double> &v, const unsigned m, const unsigned n,
                     const std::vector<double> &b, std::vector<double> &x);

} // end namespace rtt_linear

//------------------------------------------------------------------------------------------------//
// end of svbksb_pt.cc
//------------------------------------------------------------------------------------------------//
