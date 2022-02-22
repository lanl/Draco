//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   c4/scatterv.hh
 * \author Thomas M. Evans
 * \date   Thu Mar 21 11:42:03 2002
 * \brief  Data scatterv functions
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved.
 *
 * This file contains the declarations for determinate and indeterminate variable length scatter
 * functions.
 */
//------------------------------------------------------------------------------------------------//

#ifndef c4_scatterv_hh
#define c4_scatterv_hh

#include "C4_Traits.hh"
#include "C4_sys_times.h" // defines the struct DRACO_TIME_TYPE (tms for Linux).
#include <vector>

namespace rtt_c4 {
//------------------------------------------------------------------------------------------------//
/*!
 * \brief Scatter messages of known but processor-dependent size
 *
 * This subroutine handles the case where the lengths of each processor's message are known in
 * advance.
 *
 * \param outgoing_data Data to be sent from root processor. Ignored on any processor but the root
 *          processor.
 * \param incoming_data On entry, the size of each subarray must be set to the expected size of the
 *          incoming message. On return, contains the scattered data.
 */
template <typename T>
void determinate_scatterv(std::vector<std::vector<T>> &outgoing_data,
                          std::vector<T> &incoming_data);

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Scatter messages of unknown size
 *
 * This subroutine handles the case where the lengths of each processor's message are not known in
 * advance.
 *
 * \param outgoing_data Data to be sent from root processor. Ignored on all other processors.
 * \param incoming_data On return, contains the scattered data.
 */
template <typename T>
void indeterminate_scatterv(std::vector<std::vector<T>> &outgoing_data,
                            std::vector<T> &incoming_data);

} // end namespace rtt_c4

#endif // c4_scatterv_hh

//------------------------------------------------------------------------------------------------//
// end of c4/scatterv.hh
//------------------------------------------------------------------------------------------------//
