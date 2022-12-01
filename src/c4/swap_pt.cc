//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   c4/swap_pt.cc
 * \author Thomas M. Evans
 * \date   Mon Mar 25 14:44:54 2002
 * \brief  C4 MPI determinate and indeterminate swap instantiations.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "C4_Functions.hh"
#include "C4_Req.hh"
#include "swap.t.hh"
#include <c4/config.h>

namespace rtt_c4 {

//------------------------------------------------------------------------------------------------//
// EXPLICIT INSTANTIATIONS OF NON-BLOCKING SEND/RECEIVE
//------------------------------------------------------------------------------------------------//

// Unsigned

template void determinate_swap(std::vector<unsigned> const &outgoing_pid,
                               std::vector<std::vector<unsigned>> const &outgoing_data,
                               std::vector<unsigned> const &incoming_pid,
                               std::vector<std::vector<unsigned>> &incoming_data, int tag);

template void semideterminate_swap(std::vector<unsigned> const &outgoing_pid,
                                   std::vector<std::vector<unsigned>> const &outgoing_data,
                                   std::vector<unsigned> const &incoming_pid,
                                   std::vector<std::vector<unsigned>> &incoming_data, int tag);

template void determinate_swap(std::vector<std::vector<unsigned>> const &outgoing_data,
                               std::vector<std::vector<unsigned>> &incoming_data, int tag);

// Double

template void determinate_swap(std::vector<std::vector<double>> const &outgoing_data,
                               std::vector<std::vector<double>> &incoming_data, int tag);

template void determinate_swap(std::vector<unsigned> const &outgoing_pid,
                               std::vector<std::vector<double>> const &outgoing_data,
                               std::vector<unsigned> const &incoming_pid,
                               std::vector<std::vector<double>> &incoming_data, int tag);

} // end namespace rtt_c4

//------------------------------------------------------------------------------------------------//
// end of swap_pt.cc
//------------------------------------------------------------------------------------------------//
