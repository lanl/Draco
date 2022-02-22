//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   c4/swap.hh
 * \author Thomas M. Evans
 * \date   Thu Mar 21 11:42:03 2002
 * \brief  Data swap functions
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved.
 *
 * This file contains the declarations for determinate and indeterminate data swap functions.
 */
//------------------------------------------------------------------------------------------------//

#ifndef c4_swap_hh
#define c4_swap_hh

#include "C4_Traits.hh"
#include "C4_sys_times.h"
#include <vector>

namespace rtt_c4 {
//------------------------------------------------------------------------------------------------//
/*!
 * \brief Send and receive a known pattern of messages of known size.
 *
 * It is common for a program to reach a point where processes wish to exchange data with other
 * processes in a way which is hard to synchronize. This function exchanges the data using
 * asynchronous communications in the case where each process knows which other processes it is
 * expecting data from, and how much data to expect.
 *
 * \param outgoing_pid Processor ids to which this processor wishes to send data.
 *
 * \param outgoing_data Data to be send to other processors.
 *
 * \param incoming_pid Processors ids from which this processor expects to receive data.
 *
 * \param incoming_data On entry, the size of each subarray must be set to the expected size of the
 *        incoming message. On return, contains the received data.
 *
 * \param tag Tag for this exchange of data.
 */
template <typename T>
void determinate_swap(std::vector<unsigned> const &outgoing_pid,
                      std::vector<std::vector<T>> const &outgoing_data,
                      std::vector<unsigned> const &incoming_pid,
                      std::vector<std::vector<T>> &incoming_data, int tag = C4_Traits<T *>::tag);

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Send and receive a known pattern of messages of known size.
 *
 * It is common for a program to reach a point where processes wish to exchange data with other
 * processes in a way which is hard to synchronize. This function exchanges the data using
 * asynchronous communications in the case where each process knows which other processes it is
 * expecting data from, and how much data to expect.
 *
 * \param outgoing_data Data to be send to other processors. If
 *        outgoing_data[p].size()==0, then no message is sent to processor p.
 *
 * \param incoming_data On entry, the size of each subarray must be set to the
 *        expected size of the incoming message. On return, contains the
 *        received data. If incoming_data[p].size()==0, then no message is
 *        looked for from processor p.
 *
 * \param tag Tag for this exchange of data.
 */
template <typename T>
void determinate_swap(std::vector<std::vector<T>> const &outgoing_data,
                      std::vector<std::vector<T>> &incoming_data, int tag = C4_Traits<T *>::tag);

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Send and receive a known pattern of messages of unknown size.
 *
 * It is common for a program to reach a point where processes wish to exchange data with other
 * processes in a way which is hard to synchronize. This function exchanges the data using
 * asynchronous communications in the case where each process knows which other processes it is
 * expecting data from, but not how much data to expect.
 *
 * \param outgoing_pid Processor ids to which this processor wishes to send data.
 *
 * \param outgoing_data Data to be send to other processors.
 *
 * \param incoming_pid Processors ids from which this processor expects to receive data.
 *
 * \param incoming_data On return, contains the received data.
 *
 * \param tag Tag for this exchange of data.
 */
template <typename T>
void semideterminate_swap(std::vector<unsigned> const &outgoing_pid,
                          std::vector<std::vector<T>> const &outgoing_data,
                          std::vector<unsigned> const &incoming_pid,
                          std::vector<std::vector<T>> &incoming_data,
                          int tag = C4_Traits<T *>::tag);

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Send and receive an unknown pattern of messages of unknown size.
 *
 * It is common for a program to reach a point where processes wish to exchange data with other
 * processes in a way which is hard to synchronize. This function exchanges the data using
 * asynchronous communications in the case where processes do not know in advance which other
 * processes will be sending them messages..
 *
 * \param outgoing_pid Processor ids to which this processor wishes to send data.
 *
 * \param outgoing_data Data to be send to other processors.
 *
 * \param incoming_pid On return, contains processors ids from which this processor received data.
 *
 * \param incoming_data On return, contains the received data.
 *
 * \param tag Tag for this exchange of data.
 */
template <typename T>
void indeterminate_swap(std::vector<unsigned> const &outgoing_pid,
                        std::vector<std::vector<T>> const &outgoing_data,
                        std::vector<unsigned> &incoming_pid,
                        std::vector<std::vector<T>> &incoming_data, int tag = C4_Traits<T *>::tag);

} // end namespace rtt_c4

#endif // c4_swap_hh

//------------------------------------------------------------------------------------------------//
// end of c4/swap.hh
//------------------------------------------------------------------------------------------------//
