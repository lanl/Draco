//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   linear/rsolv_pt.cc
 * \author Kent Budge
 * \brief  Explicit instantiation of rsolv
 * \note   Copyright (C) 2010-2026 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "rsolv.hh" // NOLINT(misc-include-cleaner)
// #include <vector>

namespace rtt_linear {

//------------------------------------------------------------------------------------------------//
// T=vector<double> const &
//------------------------------------------------------------------------------------------------//

template void rsolv(const std::vector<double> &R, const std::size_t n, std::vector<double> &b);

} // end namespace rtt_linear

//------------------------------------------------------------------------------------------------//
// end of rsolv.cc
//------------------------------------------------------------------------------------------------//
