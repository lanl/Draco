//-----------------------------------*-C++-*----------------------------------//
/*!
* \file   c4/C4_Serial.t.hh
* \author Thomas M. Evans, Kelly Thompson <kgt@lanl.gov>
* \date   Mon Mar 25 17:06:25 2002
* \brief  Implementation of C4 serial option.
* \note   Copyright (C) 2018-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_c4_serial_t_hh
#define rtt_c4_serial_t_hh

#include "c4/config.h"

#ifdef C4_SCALAR

#include "C4_Functions.hh"
#include "C4_sys_times.h"
#include "ds++/SystemCall.hh"
#include <chrono>
#include <cstdlib>
#include <ctime>

namespace rtt_c4 {

//------------------------------------------------------------------------------------------------//
// BLOCKING SEND/RECEIVE OPERATIONS
//------------------------------------------------------------------------------------------------//

template <typename T>
int send(const T * /*buffer*/, int /*size*/, int /*destination*/, int /*tag*/) {
  return C4_SUCCESS;
}

//------------------------------------------------------------------------------------------------//
template <typename T> int receive(T * /*buffer*/, int /*size*/, int /*source*/, int /*tag*/) {
  return C4_SUCCESS;
}

//------------------------------------------------------------------------------------------------//
template <typename T>
int send_udt(const T * /*buffer*/, int /*size*/, int /*destination*/, C4_Datatype & /*data_type*/,
             int /*tag*/) {
  return C4_SUCCESS;
}

//------------------------------------------------------------------------------------------------//
template <typename T>
int receive_udt(T * /*buffer*/, int /*size*/, int /*destination*/, C4_Datatype & /*data_type*/,
                int /*tag*/) {
  return C4_SUCCESS;
}

//------------------------------------------------------------------------------------------------//
template <typename TS, typename TR>
int send_receive(TS * /*sendbuf*/, int /*sendcount*/, int /*destination*/, TR * /*recvbuf*/,
                 int /*recvcount*/, int /*source*/, int /*sendtag*/, int /*recvtag*/) {
  Insist(false, "send_receive is not support for C4_SCALAR builds.");
  return 1;
}

//------------------------------------------------------------------------------------------------//
// NON-BLOCKING SEND/RECEIVE OPERATIONS
//------------------------------------------------------------------------------------------------//

template <typename T>
C4_Req send_async(T const * /*buffer*/, int /*size*/, int /*destination*/, int /*tag*/) {
  // make a c4 request handle
  C4_Req request;
  return request;
}

//------------------------------------------------------------------------------------------------//
template <typename T>
void send_async(C4_Req &Remember(request), T const * /*buffer*/, int /*size*/, int /*destination*/,
                int /*tag*/) {
  Require(!request.inuse());
}

//------------------------------------------------------------------------------------------------//
template <typename T>
void send_is(C4_Req & /*request*/, const T * /*buffer*/, int /*size*/, int /*destination*/,
             int /*tag*/) {
  Insist(false, "send_is is not support for C4_SCALAR builds.");
}

//------------------------------------------------------------------------------------------------//
template <typename T>
C4_Req receive_async(T * /*buffer*/, int /*size*/, int /*source*/, int /*tag*/) {
  // make a c4 request handle
  C4_Req request;
  return request;
}

//------------------------------------------------------------------------------------------------//
template <typename T>
void receive_async(C4_Req &Remember(request), T * /*buffer*/, int /*size*/, int /*source*/,
                   int /*tag*/) {
  Require(!request.inuse());
}

//------------------------------------------------------------------------------------------------//
// GATHER/SCATTER
//------------------------------------------------------------------------------------------------//

template <typename T> int gather(T *send_buffer, T *receive_buffer, int size) {
  std::copy(send_buffer, send_buffer + size, receive_buffer);
  return C4_SUCCESS;
}

//------------------------------------------------------------------------------------------------//
template <typename T> int allgather(T *send_buffer, T *receive_buffer, int size) {
  std::copy(send_buffer, send_buffer + size, receive_buffer);
  return C4_SUCCESS;
}

//------------------------------------------------------------------------------------------------//
template <typename T> int scatter(T *send_buffer, T *receive_buffer, int size) {
  std::copy(send_buffer, send_buffer + size, receive_buffer);
  return C4_SUCCESS;
}

//------------------------------------------------------------------------------------------------//
template <typename T>
int gatherv(T *send_buffer, int send_size, T *receive_buffer, int * /*receive_sizes*/,
            int *receive_displs) {

  std::copy(send_buffer, send_buffer + send_size, receive_buffer + receive_displs[0]);
  return C4_SUCCESS;
}

//------------------------------------------------------------------------------------------------//
template <typename T>
int allgatherv(T *send_buffer, int send_size, T *receive_buffer, int * /*receive_sizes*/,
               int *receive_displs) {

  std::copy(send_buffer, send_buffer + send_size, receive_buffer + receive_displs[0]);
  return C4_SUCCESS;
}

//------------------------------------------------------------------------------------------------//
template <typename T>
int scatterv(T *send_buffer, int *send_sizes, int *send_displs, T *receive_buffer,
             int /*receive_size*/) {
  std::copy(send_buffer + send_displs[0], send_buffer + send_displs[0] + send_sizes[0],
            receive_buffer);

  return C4_SUCCESS;
}

//------------------------------------------------------------------------------------------------//
// GLOBAL REDUCTIONS
//------------------------------------------------------------------------------------------------//

template <typename T> void global_sum(T & /*x*/) { /* empty */
}

//------------------------------------------------------------------------------------------------//
template <typename T> void global_isum(T &send_buffer, T &receive_buffer, C4_Req & /*request*/) {
  receive_buffer = send_buffer;
}

//-----------------------------------------------------------------------------//
template <typename T> void global_prod(T * /*x*/, int /*n*/) { /* empty */
}

//------------------------------------------------------------------------------------------------//
template <typename T> void global_min(T * /*x*/, int /*n*/) { /* empty */
}

//------------------------------------------------------------------------------------------------//
template <typename T> void global_max(T * /*x*/, int /*n*/) { /* empty */
}

//------------------------------------------------------------------------------------------------//
template <typename T> void global_prod(T & /*x*/) { /* empty */
}

//------------------------------------------------------------------------------------------------//
template <typename T> void global_min(T & /*x*/) { /* empty */
}

//------------------------------------------------------------------------------------------------//
template <typename T> void global_max(T & /*x*/) { /* empty */
}

//------------------------------------------------------------------------------------------------//
template <typename T> void global_and(T & /*x*/) { /* empty */
}

//------------------------------------------------------------------------------------------------//
template <typename T> T prefix_sum(T const node_value) { return node_value; }

//------------------------------------------------------------------------------------------------//
template <typename T> void prefix_sum(T * /*buffer*/, const int32_t /*n*/) { /* empty */
}

} // end namespace rtt_c4

#endif // C4_SCALAR

#endif // rtt_c4_serial_t_hh

//------------------------------------------------------------------------------------------------//
// end of C4_Serial.t.hh
//------------------------------------------------------------------------------------------------//
