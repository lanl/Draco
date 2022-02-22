//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   c4/gatherv.hh
 * \author Thomas M. Evans
 * \date   Thu Mar 21 11:42:03 2002
 * \brief  Data gatherv functions
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved.
 *
 * This file contains the declarations for determinate and indeterminate variable length gather
 * functions.
 */
//------------------------------------------------------------------------------------------------//

#ifndef c4_gatherv_hh
#define c4_gatherv_hh

#include "C4_Traits.hh"
#include "C4_sys_times.h"
#include <vector>

namespace rtt_c4 {
//------------------------------------------------------------------------------------------------//
/*!
 * \brief Gather messages of known but processor-dependent size
 *
 * When MPI is enabled, this subroutine wraps a call to rtt_c4::gatherv to gather data across
 * processors to rank 0; it permits setting the input as a vector which can be 0 size, and
 * de-serializes the receive buffer to a vector, indexed by rank, of vectors of data.
 * "determinate" implies this form of the subroutine handles the case where the lengths of each
 *  processor's message are known in advance.
 *
 * \param[in] outgoing_data Data to be send to root processor.
 *
 * \param[in,out] incoming_data Ignored on any processor but the root processor. On the root
 *          processor, the size of each subarray must be set to the expected size of the incoming
 *          message. On return, contains the gathered data.
 */
template <class T>
void determinate_gatherv(std::vector<T> &outgoing_data, std::vector<std::vector<T>> &incoming_data);

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Gather messages of known but processor-dependent size, from all ranks to all ranks
 *
 * This subroutine performs the same operation as determinate_gatherv, but gathers the data to each
 * rank instead of just rank 0.
 *
 * \param[in] outgoing_data Data to be send to root processor.
 *
 * \param[in,out] incoming_data Gathered data per rank (outer index), available for each rank
 */
template <class T>
void determinate_allgatherv(std::vector<T> &outgoing_data,
                            std::vector<std::vector<T>> &incoming_data);

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Gather messages of unknown size
 *
 * When MPI is enabled, this subroutine wraps a call to rtt_c4::gatherv to gather data across
 * processors to rank 0; it permits setting the input as a vector which can be 0 size, and
 * de-serializes the receive buffer to a vector, indexed by rank, of vectors of data.
 * "indeterminate" implies this form of the subroutine handles the case where the lengths of each
 * processor's message are not known in advance. Consequently, an inital gather of counts and
 * displacements for the receiving buffer is needed.
 *
 * \param[in] outgoing_data Data to be send to root processor.
 *
 * \param[in,out] incoming_data Ignored on any processor but the root processor. On the root
 *           processor, on return, contains the gathered data.
 */
template <class T>
void indeterminate_gatherv(std::vector<T> &outgoing_data,
                           std::vector<std::vector<T>> &incoming_data);

void indeterminate_gatherv(std::string &outgoing_data, std::vector<std::string> &incoming_data);

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Gather messages of unknown but processor-dependent size, from all ranks to all ranks
 *
 * This subroutine performs the same operation as indeterminate_gatherv, but gathers the data to
 * each rank instead of just rank 0.
 *
 * \param[in] outgoing_data Data to be send to root processor.
 *
 * \param[in,out] incoming_data Gathered data per rank (outer index), available for each rank
 */
template <class T>
void indeterminate_allgatherv(std::vector<T> &outgoing_data,
                              std::vector<std::vector<T>> &incoming_data);

} // end namespace rtt_c4

#endif // c4_gatherv_hh

//------------------------------------------------------------------------------------------------//
// end of c4/gatherv.hh
//------------------------------------------------------------------------------------------------//
