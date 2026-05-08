//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   linear/tred2_pt.cc
 * \author Kent Budge
 * \brief  Specializations of tred2
 * \note   Copyright (C) 2010-2026 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "tred2.t.hh" // NOLINT(misc-include-cleaner)

namespace rtt_linear {

//------------------------------------------------------------------------------------------------//
// T1=T2=T3=vector<double>
//------------------------------------------------------------------------------------------------//
template void tred2(std::vector<double> &a, unsigned n, std::vector<double> &d,
                    std::vector<double> &e);

} // end namespace rtt_linear

//------------------------------------------------------------------------------------------------//
// end of linear/tred2_pt.cc
//------------------------------------------------------------------------------------------------//
