//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   linear/btridag_pt.cc
 * \author Kent Budge
 * \brief  Implementation of block tridiagonal solver
 * \note   Copyright (C) 2010-2026 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "btridag.t.hh" // NOLINT

namespace rtt_linear {

template void btridag(std::vector<double> const &a, std::vector<double> const &b,
                      std::vector<double> const &c, std::vector<double> const &r, unsigned const n,
                      unsigned const m, std::vector<double> &u);

} // end namespace rtt_linear

//------------------------------------------------------------------------------------------------//
// end of btridag_pt.cc
//------------------------------------------------------------------------------------------------//
