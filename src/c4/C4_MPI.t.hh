//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   c4/C4_MPI.t.hh
 * \author Thomas M. Evans
 * \date   Thu Mar 21 16:56:17 2002
 * \brief  C4 MPI template implementation for explicit instantiation.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef c4_C4_MPI_t_hh
#define c4_C4_MPI_t_hh

#include "C4_Functions.hh"
#include <vector>

#ifdef C4_MPI

#include "C4_Req.hh"
#include "MPI_Traits.hh"

namespace rtt_c4 {

//------------------------------------------------------------------------------------------------//
// BLOCKING SEND/RECEIVE OPERATIONS
//------------------------------------------------------------------------------------------------//

template <typename T> int send(const T *buffer, int size, int destination, int tag) {
  Require(size == 0 || buffer != nullptr);

  MPI_Send(const_cast<T *>(buffer), size, MPI_Traits<T>::element_type(), destination, tag,
           communicator);
  return C4_SUCCESS;
}

//------------------------------------------------------------------------------------------------//
template <typename T> int receive(T *buffer, int size, int source, int tag) {
  Require(size == 0 || buffer != nullptr);

  // get a handle to the MPI_Status
  MPI_Status status;

  // do the blocking receive
  Remember(int check =)
      MPI_Recv(buffer, size, MPI_Traits<T>::element_type(), source, tag, communicator, &status);
  Check(check == MPI_SUCCESS);

  // get the count of received data
  int count = 0;
  MPI_Get_count(&status, MPI_Traits<T>::element_type(), &count);
  return count;
}

//------------------------------------------------------------------------------------------------//
template <typename T>
int send_udt(const T *buffer, int size, int destination, C4_Datatype &data_type, int tag) {
  Require(size == 0 || buffer != nullptr);

  MPI_Send(const_cast<T *>(buffer), size, data_type, destination, tag, communicator);
  return C4_SUCCESS;
}

//------------------------------------------------------------------------------------------------//
template <typename T>
int receive_udt(T *buffer, int size, int source, C4_Datatype &data_type, int tag) {
  Require(size == 0 || buffer != nullptr);

  // get a handle to the MPI_Status
  MPI_Status status;

  // do the blocking receive
  MPI_Recv(buffer, size, data_type, source, tag, communicator, &status);

  // get the count of received data
  int count = 0;
  MPI_Get_count(&status, data_type, &count);
  return count;
}

//------------------------------------------------------------------------------------------------//
template <typename TS, typename TR>
int send_receive(TS *sendbuf, int sendcount, int destination, TR *recvbuf, int recvcount,
                 int source, int sendtag, int recvtag) {
  Require(sendcount == 0 || sendbuf != nullptr);
  Require(recvcount == 0 || recvbuf != nullptr);
  // buffers must not overlap
  Require(recvbuf + recvcount <= sendbuf || recvbuf >= sendbuf + sendcount);

  int check = MPI_Sendrecv(sendbuf, sendcount, MPI_Traits<TS>::element_type(), destination, sendtag,
                           recvbuf, recvcount, MPI_Traits<TR>::element_type(), source, recvtag,
                           communicator, MPI_STATUS_IGNORE);
  return check;
}

//------------------------------------------------------------------------------------------------//
// NON-BLOCKING SEND/RECEIVE OPERATIONS
//------------------------------------------------------------------------------------------------//

template <typename T> C4_Req send_async(const T *buffer, int size, int destination, int tag) {
  Require(size == 0 || buffer != nullptr);

  // make a c4 request handle
  C4_Req request;

  // do an MPI_Isend (non-blocking send)
  Remember(int const retval =)
      MPI_Isend(const_cast<T *>(buffer), size, MPI_Traits<T>::element_type(), destination, tag,
                communicator, &request.r());
  Check(retval == MPI_SUCCESS);

  // set the request to active
  request.set();
  return request;
}

//------------------------------------------------------------------------------------------------//
template <typename T>
void send_async(C4_Req &request, const T *buffer, int size, int destination, int tag) {
  Require(!request.inuse());
  Require(size == 0 || buffer != nullptr);

  // set the request
  request.set();

  // post an MPI_Isend
  MPI_Isend(const_cast<T *>(buffer), size, MPI_Traits<T>::element_type(), destination, tag,
            communicator, &request.r());
}

//------------------------------------------------------------------------------------------------//
template <typename T>
void send_is(C4_Req &request, const T *buffer, int size, int destination, int tag) {
  Require(!request.inuse());
  Require(size == 0 || buffer != nullptr);

  // set the request
  request.set();

  Remember(int const retval =)
      MPI_Issend(const_cast<T *>(buffer), size, MPI_Traits<T>::element_type(), destination, tag,
                 communicator, &request.r());
  Check(retval == MPI_SUCCESS);

  return;
}

//------------------------------------------------------------------------------------------------//
template <typename T> C4_Req receive_async(T *buffer, int size, int source, int tag) {

  Require(size == 0 || buffer != nullptr);

  // make a c4 request handle
  C4_Req request;

  // post an MPI_Irecv (non-blocking receive)
  Remember(int const retval =) MPI_Irecv(buffer, size, MPI_Traits<T>::element_type(), source, tag,
                                         communicator, &request.r());
  Check(retval == MPI_SUCCESS);

  // set the request to active
  request.set();
  return request;
}

//------------------------------------------------------------------------------------------------//
template <typename T>
void receive_async(C4_Req &request, T *buffer, int size, int source, int tag) {
  Require(size == 0 || buffer != nullptr);
  Require(!request.inuse());

  // set the request
  request.set();

  // post an MPI_Irecv
  Remember(int const retval =) MPI_Irecv(buffer, size, MPI_Traits<T>::element_type(), source, tag,
                                         communicator, &request.r());
  Check(retval == MPI_SUCCESS);
  return;
}

//------------------------------------------------------------------------------------------------//
// BROADCAST
//------------------------------------------------------------------------------------------------//

template <typename T> int broadcast(T *buffer, int size, int root) {
  Require(root >= 0 && root < nodes());
  int r = MPI_Bcast(buffer, size, MPI_Traits<T>::element_type(), root, communicator);
  return r;
}
template <typename T> int broadcast(T *buffer, size_t size, size_t root) {
  Require(root < nranks());
  Check(size < INT32_MAX);
  Check(root < INT32_MAX);
  int r = MPI_Bcast(buffer, static_cast<int>(size), MPI_Traits<T>::element_type(),
                    static_cast<int>(root), communicator);
  return r;
}

//------------------------------------------------------------------------------------------------//
// GATHER/SCATTER
//------------------------------------------------------------------------------------------------//

template <typename T> int gather(T *send_buffer, T *receive_buffer, int size) {
  int Result = MPI_Gather(send_buffer, size, MPI_Traits<T>::element_type(), receive_buffer, size,
                          MPI_Traits<T>::element_type(),
                          0, // root is always processor 0 at present
                          communicator);

  return Result;
}

//------------------------------------------------------------------------------------------------//
template <typename T> int allgather(T *send_buffer, T *receive_buffer, int size) {
  int Result = MPI_Allgather(send_buffer, size, MPI_Traits<T>::element_type(), receive_buffer, size,
                             MPI_Traits<T>::element_type(), communicator);

  return Result;
}

//------------------------------------------------------------------------------------------------//
template <typename T> int scatter(T *send_buffer, T *receive_buffer, int size) {
  int Result = MPI_Scatter(send_buffer, size, MPI_Traits<T>::element_type(), receive_buffer, size,
                           MPI_Traits<T>::element_type(),
                           0, // root is always processor 0 at present
                           communicator);

  return Result;
}

//------------------------------------------------------------------------------------------------//
template <typename T>
int gatherv(T *send_buffer, int send_size, T *receive_buffer, int *receive_sizes,
            int *receive_displs) {
  int Result = MPI_Gatherv(send_buffer, send_size, MPI_Traits<T>::element_type(), receive_buffer,
                           receive_sizes, receive_displs, MPI_Traits<T>::element_type(),
                           0, // root is always processor 0 at present
                           communicator);

  return Result;
}

//------------------------------------------------------------------------------------------------//
template <typename T>
int allgatherv(T *send_buffer, int send_size, T *receive_buffer, int *receive_sizes,
               int *receive_displs) {
  int Result =
      MPI_Allgatherv(send_buffer, send_size, MPI_Traits<T>::element_type(), receive_buffer,
                     receive_sizes, receive_displs, MPI_Traits<T>::element_type(), communicator);

  return Result;
}

//------------------------------------------------------------------------------------------------//
template <typename T>
int scatterv(T *send_buffer, int *send_sizes, int *send_displs, T *receive_buffer,
             int receive_size) {
  int Result =
      MPI_Scatterv(send_buffer, send_sizes, send_displs, MPI_Traits<T>::element_type(),
                   receive_buffer, receive_size, MPI_Traits<T>::element_type(), 0, communicator);

  return Result;
}

//------------------------------------------------------------------------------------------------//
// GLOBAL REDUCTIONS
//------------------------------------------------------------------------------------------------//

template <typename T> void global_sum(T &x) {

  // do global MPI reduction (result is on all processors) into x
  Remember(int check =)
      MPI_Allreduce(MPI_IN_PLACE, &x, 1, MPI_Traits<T>::element_type(), MPI_SUM, communicator);
  Check(check == MPI_SUCCESS);
  return;
}

//------------------------------------------------------------------------------------------------//
template <typename T> void global_isum(T &send_buffer, T &recv_buffer, C4_Req &request) {

  // do global MPI non-blocking reduction (result is on all processors) into recv_buffer
  Remember(int check =) MPI_Iallreduce(&send_buffer, &recv_buffer, 1, MPI_Traits<T>::element_type(),
                                       MPI_SUM, communicator, &(request.r()));
  request.set();
  Check(check == MPI_SUCCESS);
  return;
}

//------------------------------------------------------------------------------------------------//
template <typename T> void global_prod(T &x) {

  // do global MPI reduction (result is on all processors) into x
  MPI_Allreduce(MPI_IN_PLACE, &x, 1, MPI_Traits<T>::element_type(), MPI_PROD, communicator);
}

//------------------------------------------------------------------------------------------------//

template <typename T> void global_min(T &x) {

  // do global MPI reduction (result is on all processors) into x
  MPI_Allreduce(MPI_IN_PLACE, &x, 1, MPI_Traits<T>::element_type(), MPI_MIN, communicator);
}

//------------------------------------------------------------------------------------------------//

template <typename T> void global_max(T &x) {

  // do global MPI reduction (result is on all processors) into x
  MPI_Allreduce(MPI_IN_PLACE, &x, 1, MPI_Traits<T>::element_type(), MPI_MAX, communicator);
}

//------------------------------------------------------------------------------------------------//

template <typename T> void global_and(T &x) {

  // do global MPI reduction (result is on all processors) into x
  MPI_Allreduce(MPI_IN_PLACE, &x, 1, MPI_Traits<T>::element_type(), MPI_LAND, communicator);
}

//------------------------------------------------------------------------------------------------//

template <typename T> void global_prod(T *x, int n) {
  Require(x != nullptr);
  Require(n > 0);

  // do a element-wise global reduction (result is on all processors) into
  // x
  MPI_Allreduce(MPI_IN_PLACE, x, n, MPI_Traits<T>::element_type(), MPI_PROD, communicator);
}

//------------------------------------------------------------------------------------------------//

template <typename T> void global_min(T *x, int n) {
  Require(x != nullptr);
  Require(n > 0);

  // do a element-wise global reduction (result is on all processors) into x
  MPI_Allreduce(MPI_IN_PLACE, x, n, MPI_Traits<T>::element_type(), MPI_MIN, communicator);
}

//------------------------------------------------------------------------------------------------//

template <typename T> void global_max(T *x, int n) {
  Require(x != nullptr);
  Require(n > 0);

  // do a element-wise global reduction (result is on all processors) into x
  MPI_Allreduce(MPI_IN_PLACE, x, n, MPI_Traits<T>::element_type(), MPI_MAX, communicator);
}

//------------------------------------------------------------------------------------------------//

template <typename T> T prefix_sum(const T node_value) {
  T prefix_sum = 0;
  MPI_Scan(&node_value, &prefix_sum, 1, MPI_Traits<T>::element_type(), MPI_SUM, communicator);
  return prefix_sum;
}

//------------------------------------------------------------------------------------------------//

template <typename T> void prefix_sum(T *buffer, const int32_t n) {
  Require(buffer != nullptr);
  Require(n > 0);
  MPI_Scan(MPI_IN_PLACE, buffer, n, MPI_Traits<T>::element_type(), MPI_SUM, communicator);
}

} // end namespace rtt_c4

#endif // C4_MPI

#endif // c4_C4_MPI_t_hh

//------------------------------------------------------------------------------------------------//
// end of c4/C4_MPI.t.hh
//------------------------------------------------------------------------------------------------//
