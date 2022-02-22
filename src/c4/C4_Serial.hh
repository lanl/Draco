//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   c4/C4_Serial.hh
 * \author Thomas M. Evans
 * \date   Mon Mar 25 17:06:25 2002
 * \brief  Serial implementation of C4.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_c4_C4_Serial_hh
#define rtt_c4_C4_Serial_hh

#include "c4/config.h"
#include "ds++/Assert.hh"
#include <algorithm>

#if !defined(rtt_c4_global_hh) && !defined(rtt_c4_C4_Functions_hh)
#error "Include c4/C4_Functions.hh instead of this c4/C4_Serial.hh"
#endif

#ifdef C4_SCALAR

#include "C4_Functions.hh"
#include "C4_Req.hh"
#include "C4_Tags.hh"

namespace rtt_c4 {

//------------------------------------------------------------------------------------------------//
// SETUP FUNCTIONS
//------------------------------------------------------------------------------------------------//

template <typename Comm> void inherit(const Comm & /*comm*/) {}

template <typename T>
int create_vector_type(unsigned /*count*/, unsigned /*blocklength*/, unsigned /*stride*/,
                       C4_Datatype & /*new_type*/) {
  return C4_SUCCESS;
}

//------------------------------------------------------------------------------------------------//
// Global_<Op> functions
//------------------------------------------------------------------------------------------------//

template <typename T, typename L, typename std::enable_if<std::is_integral<L>::value, bool>::type>
void global_sum(T * /*x*/, L /*n*/) {
  return;
}

//------------------------------------------------------------------------------------------------//
// BLOCKING SEND/RECEIVE OPERATIONS
//------------------------------------------------------------------------------------------------//

template <typename T> int broadcast(T * /*buffer*/, int /*size*/, int /*root*/) {
  return C4_SUCCESS;
}
template <typename T> int broadcast(T * /*buffer*/, size_t /*size*/, size_t /*root*/) {
  return C4_SUCCESS;
}

template <typename ForwardIterator, typename OutputIterator>
void broadcast(ForwardIterator /*first*/, ForwardIterator /*last*/,
               OutputIterator /*result*/) { /* empty */
}

template <typename ForwardIterator, typename OutputIterator>
void broadcast(ForwardIterator /*first*/, ForwardIterator /*last*/, OutputIterator /*result*/,
               OutputIterator /*result_end*/) { /* empty */
}

//------------------------------------------------------------------------------------------------//
template <typename T>
int send(const T * /*buffer*/, int /*size*/, int /*destination*/, C4_Datatype & /*data_type*/,
         int /*tag*/) {
  return C4_SUCCESS;
}

//------------------------------------------------------------------------------------------------//
template <typename T>
int send_custom(const T * /*buffer*/, int /*size*/, int /*destination*/, int /*tag*/) {
  return C4_SUCCESS;
}

//------------------------------------------------------------------------------------------------//
template <typename T>
void send_is_custom(C4_Req & /* request */, const T * /*buffer*/, int /*size*/, int /*destination*/,
                    int /* tag = C4_Traits<T *>::tag */) {
  Insist(false, "send_is_custom is not support for C4_SCALAR builds.");
}

//------------------------------------------------------------------------------------------------//
template <typename T>
int receive(T * /*buffer*/, int /*size*/, int /*source*/, C4_Datatype & /*data_type*/,
            int /*tag*/) {
  return C4_SUCCESS;
}

//------------------------------------------------------------------------------------------------//
template <typename T>
int receive_custom(T * /*buffer*/, int size, int /*destination*/, int /*tag*/) {
  // expects a size of message returned
  return size;
}

//------------------------------------------------------------------------------------------------//
// NON-BLOCKING SEND/RECEIVE OPERATIONS
//------------------------------------------------------------------------------------------------//

template <typename T>
void receive_async_custom(C4_Req &Remember(request), T * /*buffer*/, int /*size*/,
                          int /*destination*/, int /*tag*/) {
  Require(!request.inuse());
}

//------------------------------------------------------------------------------------------------//
template <typename T> int message_size_custom(C4_Status /*status*/, const T & /*mpi_type*/) {
  int receive_count = 0;
  return receive_count;
}

} // end namespace rtt_c4

#endif // C4_SCALAR

#endif // rtt_c4_C4_Serial_hh

//------------------------------------------------------------------------------------------------//
// end of c4/C4_Serial.hh
//------------------------------------------------------------------------------------------------//
