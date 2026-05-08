//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   linear/tqli_pt.cc
 * \author Kent Budge
 * \brief  Specializations of tqli
 * \note   Copyright (C) 2010-2026 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "tqli.t.hh" // NOLINT(misc-include-cleaner)

namespace rtt_linear {

//------------------------------------------------------------------------------------------------//
// T=vector<double>
//------------------------------------------------------------------------------------------------//

template void tqli(std::vector<double> &d, std::vector<double> &e, const unsigned n,
                   std::vector<double> &z);

} // end namespace rtt_linear

//------------------------------------------------------------------------------------------------//
// end of tqli.cc
//------------------------------------------------------------------------------------------------//
