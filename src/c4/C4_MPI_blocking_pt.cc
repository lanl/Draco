//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   c4/C4_MPI_blocking_pt.cc
 * \author Thomas M. Evans
 * \date   Mon Mar 25 14:41:05 2002
 * \brief  C4 MPI Blocking Send/Recv instantiations.
 * \note   Copyright (C) 2016-2017 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <c4/config.h>

#ifdef C4_MPI

#include "C4_MPI.t.hh"

namespace rtt_c4 {

//---------------------------------------------------------------------------//
// EXPLICIT INSTANTIATIONS OF BLOCKING SEND/RECEIVE AND BROADCAST
//---------------------------------------------------------------------------//

template DLL_PUBLIC_c4 int send<char>(const char *, int, int, int);
template DLL_PUBLIC_c4 int send<unsigned char>(const unsigned char *, int, int,
                                               int);
template DLL_PUBLIC_c4 int send<short>(const short *, int, int, int);
template DLL_PUBLIC_c4 int send<unsigned short>(const unsigned short *, int,
                                                int, int);
template DLL_PUBLIC_c4 int send<int>(const int *, int, int, int);
template DLL_PUBLIC_c4 int send<unsigned int>(const unsigned int *, int, int,
                                              int);
template DLL_PUBLIC_c4 int send<long>(const long *, int, int, int);
template DLL_PUBLIC_c4 int send<long long>(const long long *, int, int, int);
template DLL_PUBLIC_c4 int send<unsigned long>(const unsigned long *, int, int,
                                               int);
template DLL_PUBLIC_c4 int send<unsigned long long>(const unsigned long long *,
                                                    int, int, int);
template DLL_PUBLIC_c4 int send<float>(const float *, int, int, int);
template DLL_PUBLIC_c4 int send<double>(const double *, int, int, int);
template DLL_PUBLIC_c4 int send<long double>(const long double *, int, int,
                                             int);

template DLL_PUBLIC_c4 int send_udt<double>(const double *, int, int,
                                            C4_Datatype &, int);

template DLL_PUBLIC_c4 int receive<char>(char *, int, int, int);
template DLL_PUBLIC_c4 int receive<unsigned char>(unsigned char *, int, int,
                                                  int);
template DLL_PUBLIC_c4 int receive<short>(short *, int, int, int);
template DLL_PUBLIC_c4 int receive<unsigned short>(unsigned short *, int, int,
                                                   int);
template DLL_PUBLIC_c4 int receive<int>(int *, int, int, int);
template DLL_PUBLIC_c4 int receive<unsigned int>(unsigned int *, int, int, int);
template DLL_PUBLIC_c4 int receive<long>(long *, int, int, int);
template DLL_PUBLIC_c4 int receive<long long>(long long *, int, int, int);
template DLL_PUBLIC_c4 int receive<unsigned long>(unsigned long *, int, int,
                                                  int);
template DLL_PUBLIC_c4 int receive<unsigned long long>(unsigned long long *,
                                                       int, int, int);
template DLL_PUBLIC_c4 int receive<float>(float *, int, int, int);
template DLL_PUBLIC_c4 int receive<double>(double *, int, int, int);
template DLL_PUBLIC_c4 int receive<long double>(long double *, int, int, int);

template DLL_PUBLIC_c4 int receive_udt<double>(double *, int, int,
                                               C4_Datatype &, int);

template DLL_PUBLIC_c4 int broadcast<char>(char *, int, int);
template DLL_PUBLIC_c4 int broadcast<unsigned char>(unsigned char *, int, int);
template DLL_PUBLIC_c4 int broadcast<short>(short *, int, int);
template DLL_PUBLIC_c4 int broadcast<unsigned short>(unsigned short *, int,
                                                     int);
template DLL_PUBLIC_c4 int broadcast<int>(int *, int, int);
template DLL_PUBLIC_c4 int broadcast<unsigned int>(unsigned int *, int, int);
template DLL_PUBLIC_c4 int broadcast<long>(long *, int, int);
template DLL_PUBLIC_c4 int broadcast<long long>(long long *, int, int);
template DLL_PUBLIC_c4 int broadcast<unsigned long>(unsigned long *, int, int);
template DLL_PUBLIC_c4 int broadcast<unsigned long long>(unsigned long long *,
                                                         int, int);
template DLL_PUBLIC_c4 int broadcast<float>(float *, int, int);
template DLL_PUBLIC_c4 int broadcast<double>(double *, int, int);
template DLL_PUBLIC_c4 int broadcast<long double>(long double *, int, int);

} // end namespace rtt_c4

#endif // C4_MPI

//---------------------------------------------------------------------------//
// end of C4_MPI_blocking_pt.cc
//---------------------------------------------------------------------------//
