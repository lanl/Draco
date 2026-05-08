//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   linear/ludcmp_pt.cc
 * \author Kent Budge
 * \brief  LU decomposition
 * \note   Copyright (C) 2010-2026 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "dsxx/Slice.hh" // NOLINT(misc-include-cleaner)
#include "ludcmp.hh"     // NOLINT(misc-include-cleaner)

namespace rtt_linear {

template void ludcmp(std::vector<double> &a, std::vector<unsigned> &indx, double &d);

template void lubksb(std::vector<double> const &a, std::vector<unsigned> const &indx,
                     std::vector<double> &b);

template void lubksb(std::vector<double> const &a, std::vector<unsigned> const &indx,
                     rtt_dsxx::Slice<std::vector<double>::iterator> &b);

} // end namespace rtt_linear

//------------------------------------------------------------------------------------------------//
// end of ludcmp_pt.cc
//------------------------------------------------------------------------------------------------//
