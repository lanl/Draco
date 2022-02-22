//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   c4/Invert_Comm_Map.hh
 * \author Mike Buksas, Rob Lowrie
 * \date   Mon Nov 19 10:09:10 2007
 * \brief  Implementation of Invert_Comm_Map
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef c4_Invert_Comm_Map_hh
#define c4_Invert_Comm_Map_hh

#include "C4_Functions.hh"
#include <map>

namespace rtt_c4 {

//! Map type for Invert_Comm_Map functions
using Invert_Comm_Map_t = std::map<int, size_t>;

//------------------------------------------------------------------------------------------------//
/**
 * \brief Invert the contents of a one-to-many mapping between nodes.
 *
 * \param[in] to_map On input, a map from processor number to the size of information to be sent to
 *        that processor by the current processor.
 * \param[out] from_map On output, a map from processor number to the size of information to be
 *        received from that processor by the current processor.  On input, ignored and deleted.
 *
 * Here, the units of the "size of information" is up to the caller.  For example, it might be the
 * number of bytes, or the number of elements in an array.  The size must be positive (specifically,
 * nonzero).
 */
void invert_comm_map(Invert_Comm_Map_t const &to_map, Invert_Comm_Map_t &from_map);

//------------------------------------------------------------------------------------------------//
/**
 * \brief Returns the number of remote processors sending to this proc.
 *
 * \param[in] first First value of input iterator over processor numbers that this processor is
 *            sending to.
 * \param[in] last Last value of input iterator.
 *
 * Note that only remote processors are counted; if this processor is sending to itself, this proc
 * is not counted as a sending proc.  Therefore, the result of this function is the number of async
 * receives that the caller should post.
 *
 * To be more useful outside of invert_comm_map(), this function could be templated on the iterator
 * type and a functor to access the value, as opposed to being restricted to the std::map iterator
 * and assuming the processor numbers are stored in iterator->first.  But invert_comm_map doesn't
 * need that, and we decided not to add that complexity.
 */
int get_num_recv(Invert_Comm_Map_t::const_iterator first, Invert_Comm_Map_t::const_iterator last);

} // end namespace rtt_c4

#endif // c4_Invert_Comm_Map_hh

//------------------------------------------------------------------------------------------------//
// end of c4/Invert_Comm_Map.hh
//------------------------------------------------------------------------------------------------//
