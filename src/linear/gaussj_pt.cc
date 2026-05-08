//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   linear/gaussj_pt.cc
 * \author Kent Budge
 * \brief  Specializations of gaussj
 * \note   Copyright (C) 2010-2026 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "gaussj.t.hh" // NOLINT(misc-include-cleaner)

namespace rtt_linear {

template void gaussj(std::vector<double> &A, unsigned n, std::vector<double> &b, unsigned m);
template void gaussj(std::vector<std::vector<double>> &A, std::vector<double> &b);

} // end namespace rtt_linear

//------------------------------------------------------------------------------------------------//
// end of linear/gaussj_pt.cc
//------------------------------------------------------------------------------------------------//
