//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   linear/btridag_pt.cc
 * \author Kent Budge
 * \date   Wed Sep 15 13:03:41 MDT 2010
 * \brief  Implementation of block tridiagonal solver
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "btridag.t.hh"
#include <vector>

namespace rtt_linear {

template void btridag(std::vector<double> const &a, std::vector<double> const &b,
                      std::vector<double> const &c, std::vector<double> const &r, unsigned const n,
                      unsigned const m, std::vector<double> &u);

} // end namespace rtt_linear

//------------------------------------------------------------------------------------------------//
// end of btridag_pt.cc
//------------------------------------------------------------------------------------------------//
