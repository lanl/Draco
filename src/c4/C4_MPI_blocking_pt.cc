//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   c4/C4_MPI_blocking_pt.cc
 * \author Thomas M. Evans
 * \date   Mon Mar 25 14:41:05 2002
 * \brief  C4 MPI Blocking Send/Receive instantiations.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "c4/config.h"

#ifdef C4_MPI
#include "C4_MPI.t.hh"
#else
#include "C4_Serial.t.hh"
#endif

namespace rtt_c4 {

//------------------------------------------------------------------------------------------------//
// EXPLICIT INSTANTIATIONS OF BLOCKING SEND/RECEIVE AND BROADCAST
//------------------------------------------------------------------------------------------------//

template int send<bool>(const bool *, int, int, int);
template int send<char>(const char *, int, int, int);
template int send<unsigned char>(const unsigned char *, int, int, int);
template int send<short>(const short *, int, int, int);
template int send<unsigned short>(const unsigned short *, int, int, int);
template int send<int>(const int *, int, int, int);
template int send<unsigned int>(const unsigned int *, int, int, int);
template int send<long>(const long *, int, int, int);
template int send<long long>(const long long *, int, int, int);
template int send<unsigned long>(const unsigned long *, int, int, int);
template int send<unsigned long long>(const unsigned long long *, int, int, int);
template int send<float>(const float *, int, int, int);
template int send<double>(const double *, int, int, int);
template int send<long double>(const long double *, int, int, int);

template int send_udt<double>(const double *, int, int, C4_Datatype &, int);

template int receive<bool>(bool *, int, int, int);
template int receive<char>(char *, int, int, int);
template int receive<unsigned char>(unsigned char *, int, int, int);
template int receive<short>(short *, int, int, int);
template int receive<unsigned short>(unsigned short *, int, int, int);
template int receive<int>(int *, int, int, int);
template int receive<unsigned int>(unsigned int *, int, int, int);
template int receive<long>(long *, int, int, int);
template int receive<long long>(long long *, int, int, int);
template int receive<unsigned long>(unsigned long *, int, int, int);
template int receive<unsigned long long>(unsigned long long *, int, int, int);
template int receive<float>(float *, int, int, int);
template int receive<double>(double *, int, int, int);
template int receive<long double>(long double *, int, int, int);

template int receive_udt<double>(double *, int, int, C4_Datatype &, int);

template int broadcast<bool>(bool *, int, int);
template int broadcast<char>(char *, int, int);
template int broadcast<unsigned char>(unsigned char *, int, int);
template int broadcast<short>(short *, int, int);
template int broadcast<unsigned short>(unsigned short *, int, int);
template int broadcast<int>(int *, int, int);
template int broadcast<unsigned int>(unsigned int *, int, int);
template int broadcast<long>(long *, int, int);
template int broadcast<long long>(long long *, int, int);
template int broadcast<unsigned long>(unsigned long *, int, int);
template int broadcast<unsigned long long>(unsigned long long *, int, int);
template int broadcast<float>(float *, int, int);
template int broadcast<double>(double *, int, int);
template int broadcast<long double>(long double *, int, int);

template int broadcast<bool>(bool *, size_t, size_t);
template int broadcast<char>(char *, size_t, size_t);
template int broadcast<unsigned char>(unsigned char *, size_t, size_t);
template int broadcast<short>(short *, size_t, size_t);
template int broadcast<unsigned short>(unsigned short *, size_t, size_t);
template int broadcast<int>(int *, size_t, size_t);
template int broadcast<unsigned int>(unsigned int *, size_t, size_t);
template int broadcast<long>(long *, size_t, size_t);
template int broadcast<long long>(long long *, size_t, size_t);
template int broadcast<unsigned long>(unsigned long *, size_t, size_t);
template int broadcast<unsigned long long>(unsigned long long *, size_t, size_t);
template int broadcast<float>(float *, size_t, size_t);
template int broadcast<double>(double *, size_t, size_t);
template int broadcast<long double>(long double *, size_t, size_t);

template int send_receive(bool *sendbuf, int sendcount, int destination, bool *recvbuf,
                          int recvcount, int source, int sendtag, int recvtag);
template int send_receive(char *sendbuf, int sendcount, int destination, char *recvbuf,
                          int recvcount, int source, int sendtag, int recvtag);
template int send_receive(int *sendbuf, int sendcount, int destination, int *recvbuf, int recvcount,
                          int source, int sendtag, int recvtag);
template int send_receive(long *sendbuf, int sendcount, int destination, long *recvbuf,
                          int recvcount, int source, int sendtag, int recvtag);
template int send_receive(float *sendbuf, int sendcount, int destination, float *recvbuf,
                          int recvcount, int source, int sendtag, int recvtag);
template int send_receive(double *sendbuf, int sendcount, int destination, double *recvbuf,
                          int recvcount, int source, int sendtag, int recvtag);

template int prefix_sum(const int node_value);
template uint32_t prefix_sum(const uint32_t node_value);
template long prefix_sum(const long node_value);
template long long prefix_sum(const long long node_value);
template uint64_t prefix_sum(const uint64_t node_value);
template float prefix_sum(const float node_value);
template double prefix_sum(const double node_value);

template void prefix_sum(int32_t *buffer, int32_t n);
template void prefix_sum(uint32_t *buffer, int32_t n);
template void prefix_sum(int64_t *buffer, int32_t n);
template void prefix_sum(uint64_t *buffer, int32_t n);
template void prefix_sum(float *buffer, int32_t n);
template void prefix_sum(double *buffer, int32_t n);

} // end namespace rtt_c4

//------------------------------------------------------------------------------------------------//
// end of C4_MPI_blocking_pt.cc
//------------------------------------------------------------------------------------------------//
