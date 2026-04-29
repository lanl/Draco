//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   c4/C4_MPI.i.hh
 * \author Alex R Long
 * \date   Mon Aug 21 07:47:01 2017
 * \brief  C4 MPI standard implementations.
 * \note   Copyright (C) 2017-2026 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef c4_C4_MPI_i_hh
#define c4_C4_MPI_i_hh

#include "C4_Req.hh"
#include <concepts>

#ifdef C4_MPI

namespace rtt_c4 {

//------------------------------------------------------------------------------------------------//
template <typename T>
void send_is_custom(C4_Req &request, const T *buffer, int size, int destination,
                    int tag /* = C4_Traits<T *>::tag */) {
  Require(!request.inuse());
  Require(buffer != nullptr);

  // set the request
  request.set();

  Remember(int const retval =) MPI_Issend(const_cast<T *>(buffer), size, T::MPI_Type, destination,
                                          tag, communicator, &request.r());
  Check(retval == MPI_SUCCESS);

  return;
}

//------------------------------------------------------------------------------------------------//
template <typename T> int send_custom(const T *buffer, int size, int destination, int tag) {
  Require(buffer != nullptr);
  MPI_Send(const_cast<T *>(buffer), size, T::MPI_Type, destination, tag, communicator);
  return C4_SUCCESS;
}

//------------------------------------------------------------------------------------------------//
template <typename T>
void receive_async_custom(C4_Req &request, T *buffer, int size, int source, int tag) {
  Require(!request.inuse());
  Require(buffer != nullptr);
  Remember(int custom_mpi_type_size);
  Remember(MPI_Type_size(T::MPI_Type, &custom_mpi_type_size));
  Require(custom_mpi_type_size == sizeof(T));

  // set the request
  request.set();

  // post an MPI_Irecv
  Remember(int const retval =)
      MPI_Irecv(buffer, size, T::MPI_Type, source, tag, communicator, &request.r());
  Check(retval == MPI_SUCCESS);
  return;
}

//------------------------------------------------------------------------------------------------//
template <typename T> int receive_custom(T *buffer, int size, int source, int tag) {
  Require(buffer != nullptr);
  // get a handle to the MPI_Status
  MPI_Status status;

  // do the blocking receive
  Remember(int check =) MPI_Recv(buffer, size, T::MPI_Type, source, tag, communicator, &status);
  Check(check == MPI_SUCCESS);

  // get the count of received data
  int count = 0;
  MPI_Get_count(&status, T::MPI_Type, &count);
  return count;
}

//------------------------------------------------------------------------------------------------//
template <typename T> int message_size_custom(C4_Status status, const T &mpi_type) {
  int receive_count = 0;
  MPI_Get_count(status.get_status_obj(), mpi_type, &receive_count);
  return receive_count;
}

} // end namespace rtt_c4

#endif // C4_MPI

#endif // c4_C4_MPI_i_hh

//------------------------------------------------------------------------------------------------//
// end of c4/C4_MPI.i.hh
//------------------------------------------------------------------------------------------------//
