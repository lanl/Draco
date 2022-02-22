//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   c4/C4_Traits.hh
 * \author Thomas M. Evans
 * \date   Thu Mar 21 16:37:29 2002
 * \brief  Traits for C4 intrinsic types.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_c4_C4_Traits_hh
#define rtt_c4_C4_Traits_hh

#include "C4_Tags.hh"

namespace rtt_c4 {

//================================================================================================//
/*!
 * \struct C4_Traits
 *
 * This struct and its specializations are used to implement the type-safe default message tags in
 * C4.  Any other type-determined property needed in C4 would also go here.
 */
//================================================================================================//

template <typename T> struct C4_Traits {};

//------------------------------------------------------------------------------------------------//
// SPECIALIZATION OF INTRINSIC ELEMENTAL TYPES
//------------------------------------------------------------------------------------------------//

template <> struct C4_Traits<bool> { static const int tag = 430; };

template <> struct C4_Traits<char> { static const int tag = 431; };

template <> struct C4_Traits<unsigned char> { static const int tag = 432; };

template <> struct C4_Traits<short> { static const int tag = 433; };

template <> struct C4_Traits<unsigned short> { static const int tag = 434; };

template <> struct C4_Traits<int> { static const int tag = 435; };

template <> struct C4_Traits<unsigned int> { static const int tag = 436; };

template <> struct C4_Traits<long> { static const int tag = 437; };

template <> struct C4_Traits<unsigned long> { static const int tag = 438; };

template <> struct C4_Traits<float> { static const int tag = 439; };

template <> struct C4_Traits<double> { static const int tag = 440; };

template <> struct C4_Traits<long double> { static const int tag = 441; };

template <> struct C4_Traits<unsigned long long> { static const int tag = 442; };

template <> struct C4_Traits<long long> { static const int tag = 443; };

//------------------------------------------------------------------------------------------------//
// SPECIALIZATION OF INTRINSIC POINTER TYPES
//------------------------------------------------------------------------------------------------//

template <> struct C4_Traits<bool *> { static const int tag = 450; };

template <> struct C4_Traits<char *> { static const int tag = 451; };

template <> struct C4_Traits<unsigned char *> { static const int tag = 452; };

template <> struct C4_Traits<short *> { static const int tag = 453; };

template <> struct C4_Traits<unsigned short *> { static const int tag = 454; };

template <> struct C4_Traits<int *> { static const int tag = 455; };

template <> struct C4_Traits<unsigned int *> { static const int tag = 456; };

template <> struct C4_Traits<long *> { static const int tag = 457; };

template <> struct C4_Traits<unsigned long *> { static const int tag = 458; };

template <> struct C4_Traits<float *> { static const int tag = 459; };

template <> struct C4_Traits<double *> { static const int tag = 460; };

template <> struct C4_Traits<long double *> { static const int tag = 461; };

template <> struct C4_Traits<unsigned long long *> { static const int tag = 462; };

template <> struct C4_Traits<long long *> { static const int tag = 463; };

} // end namespace rtt_c4

#endif // rtt_c4_C4_Traits_hh

//------------------------------------------------------------------------------------------------//
// end of c4/C4_Traits.hh
//------------------------------------------------------------------------------------------------//
