//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   c4/MPI_Traits.hh
 * \author Thomas M. Evans
 * \date   Thu Mar 21 11:07:40 2002
 * \brief  Traits classes for MPI types.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef c4_MPI_Traits_hh
#define c4_MPI_Traits_hh

#include "c4_mpi.h"

namespace rtt_c4 {

//================================================================================================//
/*!
 * \struct MPI_Traits
 *
 * \brief Provide a generic way to get MPI_Datatype arguments for MPI function
 *        calls.
 *
 * This struct provides a generic programming--common way to get MPI_Datatype arguments for MPI
 * function calls. The static function, element_type(), returns an argument of type MPI_Datatype
 * that matches a C++ datatype with an MPI_Datatype.
 */
//================================================================================================//

template <typename T> struct MPI_Traits {};

//------------------------------------------------------------------------------------------------//
// SPECIALIZATIONS OF MPI_Traits FOR DIFFERENT TYPES
//------------------------------------------------------------------------------------------------//

#ifdef C4_MPI

template <> struct MPI_Traits<bool> {
  static MPI_Datatype element_type() { return MPI_C_BOOL; }
};

template <> struct MPI_Traits<char> {
  static MPI_Datatype element_type() { return MPI_CHAR; }
};

template <> struct MPI_Traits<unsigned char> {
  static MPI_Datatype element_type() { return MPI_UNSIGNED_CHAR; }
};

template <> struct MPI_Traits<short> {
  static MPI_Datatype element_type() { return MPI_SHORT; }
};

template <> struct MPI_Traits<unsigned short> {
  static MPI_Datatype element_type() { return MPI_UNSIGNED_SHORT; }
};

template <> struct MPI_Traits<int> {
  static MPI_Datatype element_type() { return MPI_INT; }
};

template <> struct MPI_Traits<unsigned int> {
  static MPI_Datatype element_type() { return MPI_UNSIGNED; }
};

template <> struct MPI_Traits<long> {
  static MPI_Datatype element_type() { return MPI_LONG; }
};

template <> struct MPI_Traits<long long> {
  static MPI_Datatype element_type() { return MPI_LONG_LONG; }
};

template <> struct MPI_Traits<unsigned long> {
  static MPI_Datatype element_type() { return MPI_UNSIGNED_LONG; }
};

template <> struct MPI_Traits<unsigned long long> {
  static MPI_Datatype element_type() { return MPI_UNSIGNED_LONG_LONG; }
};

template <> struct MPI_Traits<float> {
  static MPI_Datatype element_type() { return MPI_FLOAT; }
};

template <> struct MPI_Traits<double> {
  static MPI_Datatype element_type() { return MPI_DOUBLE; }
};

template <> struct MPI_Traits<long double> {
  static MPI_Datatype element_type() { return MPI_LONG_DOUBLE; }
};

#endif

} // end namespace rtt_c4

#endif // c4_MPI_Traits_hh

//------------------------------------------------------------------------------------------------//
// end of c4/MPI_Traits.hh
//------------------------------------------------------------------------------------------------//
