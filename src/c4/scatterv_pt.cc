//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   c4/scatterv_pt.cc
 * \author Thomas M. Evans
 * \brief  C4 MPI determinate and indeterminate scatterv instantiations.
 * \note   Copyright (C) 2010-2026 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "scatterv.t.hh" // NOLINT(misc-include-cleaner)

namespace rtt_c4 {

//------------------------------------------------------------------------------------------------//
// EXPLICIT INSTANTIATIONS OF NON-BLOCKING SEND/RECEIVE
//------------------------------------------------------------------------------------------------//

template void indeterminate_scatterv(std::vector<std::vector<unsigned>> &outgoing_data,
                                     std::vector<unsigned> &incoming_data);
template void indeterminate_scatterv(std::vector<std::vector<int>> &outgoing_data,
                                     std::vector<int> &incoming_data);
template void indeterminate_scatterv(std::vector<std::vector<long long>> &outgoing_data,
                                     std::vector<long long> &incoming_data);
template void indeterminate_scatterv(std::vector<std::vector<double>> &outgoing_data,
                                     std::vector<double> &incoming_data);

//------------------------------------------------------------------------------------------------//
template void determinate_scatterv(std::vector<std::vector<unsigned>> &outgoing_data,
                                   std::vector<unsigned> &incoming_data);
template void determinate_scatterv(std::vector<std::vector<int>> &outgoing_data,
                                   std::vector<int> &incoming_data);
template void determinate_scatterv(std::vector<std::vector<long long>> &outgoing_data,
                                   std::vector<long long> &incoming_data);
template void determinate_scatterv(std::vector<std::vector<double>> &outgoing_data,
                                   std::vector<double> &incoming_data);

} // end namespace rtt_c4

//------------------------------------------------------------------------------------------------//
// end of C4_MPI_scatterv_pt.cc
//------------------------------------------------------------------------------------------------//
