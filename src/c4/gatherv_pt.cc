//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   c4/gatherv_pt.cc
 * \author Thomas M. Evans
 * \brief  C4 MPI determinate and indeterminate gatherv instantiations.
 * \note   Copyright (C) 2010-2026 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "gatherv.t.hh" // NOLINT(misc-include-cleaner)

namespace rtt_c4 {

//------------------------------------------------------------------------------------------------//
// EXPLICIT INSTANTIATIONS OF NON-BLOCKING SEND/RECEIVE
//------------------------------------------------------------------------------------------------//
template void indeterminate_gatherv<uint64_t>(std::vector<uint64_t> &outgoing_data,
                                              std::vector<std::vector<uint64_t>> &incoming_data,
                                              bool is_allgatherv);

template void indeterminate_gatherv<int64_t>(std::vector<int64_t> &outgoing_data,
                                             std::vector<std::vector<int64_t>> &incoming_data,
                                             bool is_allgatherv);

template void indeterminate_gatherv<unsigned>(std::vector<unsigned> &outgoing_data,
                                              std::vector<std::vector<unsigned>> &incoming_data,
                                              bool is_allgatherv);

template void indeterminate_gatherv<int>(std::vector<int> &outgoing_data,
                                         std::vector<std::vector<int>> &incoming_data,
                                         bool is_allgatherv);

template void indeterminate_gatherv<long long>(std::vector<long long> &outgoing_data,
                                               std::vector<std::vector<long long>> &incoming_data,
                                               bool is_allgatherv);

template void indeterminate_gatherv<double>(std::vector<double> &outgoing_data,
                                            std::vector<std::vector<double>> &incoming_data,
                                            bool is_allgatherv);

template void indeterminate_gatherv<char>(std::vector<char> &outgoing_data,
                                          std::vector<std::vector<char>> &incoming_data,
                                          bool is_allgatherv);

//------------------------------------------------------------------------------------------------//
template void determinate_gatherv<uint64_t>(std::vector<uint64_t> &outgoing_data,
                                            std::vector<std::vector<uint64_t>> &incoming_data,
                                            bool is_allgatherv);

template void determinate_gatherv<int64_t>(std::vector<int64_t> &outgoing_data,
                                           std::vector<std::vector<int64_t>> &incoming_data,
                                           bool is_allgatherv);

template void determinate_gatherv<unsigned>(std::vector<unsigned> &outgoing_data,
                                            std::vector<std::vector<unsigned>> &incoming_data,
                                            bool is_allgatherv);

template void determinate_gatherv<int>(std::vector<int> &outgoing_data,
                                       std::vector<std::vector<int>> &incoming_data,
                                       bool is_allgatherv);

template void determinate_gatherv<long long>(std::vector<long long> &outgoing_data,
                                             std::vector<std::vector<long long>> &incoming_data,
                                             bool is_allgatherv);

template void determinate_gatherv<double>(std::vector<double> &outgoing_data,
                                          std::vector<std::vector<double>> &incoming_data,
                                          bool is_allgatherv);

template void determinate_gatherv<char>(std::vector<char> &outgoing_data,
                                        std::vector<std::vector<char>> &incoming_data,
                                        bool is_allgatherv);

//------------------------------------------------------------------------------------------------//
template void indeterminate_allgatherv<uint64_t>(std::vector<uint64_t> &outgoing_data,
                                                 std::vector<std::vector<uint64_t>> &incoming_data);
template void indeterminate_allgatherv<int64_t>(std::vector<int64_t> &outgoing_data,
                                                std::vector<std::vector<int64_t>> &incoming_data);
template void indeterminate_allgatherv<unsigned>(std::vector<unsigned> &outgoing_data,
                                                 std::vector<std::vector<unsigned>> &incoming_data);

template void determinate_allgatherv<uint64_t>(std::vector<uint64_t> &outgoing_data,
                                               std::vector<std::vector<uint64_t>> &incoming_data);
template void determinate_allgatherv<int64_t>(std::vector<int64_t> &outgoing_data,
                                              std::vector<std::vector<int64_t>> &incoming_data);
template void determinate_allgatherv<unsigned>(std::vector<unsigned> &outgoing_data,
                                               std::vector<std::vector<unsigned>> &incoming_data);
template void determinate_allgatherv<double>(std::vector<double> &outgoing_data,
                                             std::vector<std::vector<double>> &incoming_data);

} // end namespace rtt_c4

//------------------------------------------------------------------------------------------------//
// end of C4_MPI_gatherv_pt.cc
//------------------------------------------------------------------------------------------------//
