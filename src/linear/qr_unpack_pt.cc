//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   linear/qr_unpack_pt.cc
 * \author Kent Budge
 * \brief  Specializations of qr_unpack
 * \note   Copyright (C) 2010-2026 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "qr_unpack.t.hh" // NOLINT(misc-include-cleaner)

namespace rtt_linear {

//------------------------------------------------------------------------------------------------//
// RandomContainer = vector<double>
//------------------------------------------------------------------------------------------------//

template void qr_unpack(std::vector<double> &r, const unsigned n, const std::vector<double> &c,
                        const std::vector<double> &d, std::vector<double> &qt);

} // end namespace rtt_linear

//------------------------------------------------------------------------------------------------//
// end of qr_unpack.cc
//------------------------------------------------------------------------------------------------//
