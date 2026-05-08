//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   linear/qrupdt_pt.cc
 * \author Kent Budge
 * \brief  Specializations of qrupdt
 * \note   Copyright (C) 2004-2026 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "qrupdt.t.hh" // NOLINT(misc-include-cleaner)

namespace rtt_linear {

//------------------------------------------------------------------------------------------------//
// T = RandomContainer = vector<double>
//------------------------------------------------------------------------------------------------//

template void qrupdt(std::vector<double> &r, std::vector<double> &qt, const unsigned n,
                     std::vector<double> &u, std::vector<double> &v);

} // end namespace rtt_linear

//------------------------------------------------------------------------------------------------//
//  end of linear/qrupdt_pt.cc
//------------------------------------------------------------------------------------------------//
