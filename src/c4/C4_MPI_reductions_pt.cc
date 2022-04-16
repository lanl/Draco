//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   c4/C4_MPI_reductions_pt.cc
 * \author Thomas M. Evans
 * \date   Mon Mar 25 11:12:35 2002
 * \brief  C4 MPI global reduction instantiations.
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
// EXPLICIT INSTANTIATIONS OF GLOBAL REDUCTIONS
//------------------------------------------------------------------------------------------------//

template void global_sum<short>(short &);
template void global_sum<unsigned short>(unsigned short &);
template void global_sum<int>(int &);
template void global_sum<unsigned int>(unsigned int &);
template void global_sum<long>(long &);
template void global_sum<unsigned long>(unsigned long &);
template void global_sum<float>(float &);
template void global_sum<double>(double &);
template void global_sum<long double>(long double &);
template void global_sum<long long>(long long &);
template void global_sum<unsigned long long>(unsigned long long &);

template void global_isum<short>(short &, short &, C4_Req &);
template void global_isum<unsigned short>(unsigned short &, unsigned short &, C4_Req &);
template void global_isum<int>(int &, int &, C4_Req &);
template void global_isum<unsigned int>(unsigned int &, unsigned int &, C4_Req &);
template void global_isum<long>(long &, long &, C4_Req &);
template void global_isum<unsigned long>(unsigned long &, unsigned long &, C4_Req &);
template void global_isum<float>(float &, float &, C4_Req &);
template void global_isum<double>(double &, double &, C4_Req &);
template void global_isum<long double>(long double &, long double &, C4_Req &);
template void global_isum<long long>(long long &, long long &, C4_Req &);
template void global_isum<unsigned long long>(unsigned long long &, unsigned long long &, C4_Req &);

template void global_prod<short>(short &);
template void global_prod<unsigned short>(unsigned short &);
template void global_prod<int>(int &);
template void global_prod<unsigned int>(unsigned int &);
template void global_prod<long>(long &);
template void global_prod<unsigned long>(unsigned long &);
template void global_prod<float>(float &);
template void global_prod<double>(double &);
template void global_prod<long double>(long double &);
template void global_prod<long long>(long long &);
template void global_prod<unsigned long long>(unsigned long long &);

template void global_max<short>(short &);
template void global_max<unsigned short>(unsigned short &);
template void global_max<int>(int &);
template void global_max<unsigned int>(unsigned int &);
template void global_max<long>(long &);
template void global_max<unsigned long>(unsigned long &);
template void global_max<float>(float &);
template void global_max<double>(double &);
template void global_max<long double>(long double &);
template void global_max<long long>(long long &);
template void global_max<unsigned long long>(unsigned long long &);

template void global_min<short>(short &);
template void global_min<unsigned short>(unsigned short &);
template void global_min<int>(int &);
template void global_min<unsigned int>(unsigned int &);
template void global_min<long>(long &);
template void global_min<unsigned long>(unsigned long &);
template void global_min<float>(float &);
template void global_min<double>(double &);
template void global_min<long double>(long double &);
template void global_min<long long>(long long &);
template void global_min<unsigned long long>(unsigned long long &);

template void global_and<bool>(bool &);

template void global_prod<short>(short *, int);
template void global_prod<unsigned short>(unsigned short *, int);
template void global_prod<int>(int *, int);
template void global_prod<unsigned int>(unsigned int *, int);
template void global_prod<long>(long *, int);
template void global_prod<unsigned long>(unsigned long *, int);
template void global_prod<float>(float *, int);
template void global_prod<double>(double *, int);
template void global_prod<long double>(long double *, int);
template void global_prod<long long>(long long *, int);
template void global_prod<unsigned long long>(unsigned long long *, int);

template void global_max<short>(short *, int);
template void global_max<unsigned short>(unsigned short *, int);
template void global_max<int>(int *, int);
template void global_max<unsigned int>(unsigned int *, int);
template void global_max<long>(long *, int);
template void global_max<unsigned long>(unsigned long *, int);
template void global_max<float>(float *, int);
template void global_max<double>(double *, int);
template void global_max<long double>(long double *, int);
template void global_max<long long>(long long *, int);
template void global_max<unsigned long long>(unsigned long long *, int);

template void global_min<short>(short *, int);
template void global_min<unsigned short>(unsigned short *, int);
template void global_min<int>(int *, int);
template void global_min<unsigned int>(unsigned int *, int);
template void global_min<long>(long *, int);
template void global_min<unsigned long>(unsigned long *, int);
template void global_min<float>(float *, int);
template void global_min<double>(double *, int);
template void global_min<long double>(long double *, int);
template void global_min<long long>(long long *, int);
template void global_min<unsigned long long>(unsigned long long *, int);

} // end namespace rtt_c4

//------------------------------------------------------------------------------------------------//
// end of C4_MPI_reductions_pt.cc
//------------------------------------------------------------------------------------------------//
