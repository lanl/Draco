//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   c4/MPI_Traits.hh
 * \author Thomas M. Evans
 * \date   Thu Mar 21 11:07:40 2002
 * \brief  Traits classes for MPI types.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __c4_MPI_Traits_hh__
#define __c4_MPI_Traits_hh__

#include "c4_mpi.h"

namespace rtt_c4
{
 
//===========================================================================//
/*!
 * \struct MPI_Traits
 *
 * \brief Provide a generic way to get MPI_Datatype arguments for MPI function
 * calls. 
 *
 * This struct provides a generic programming--common way to get MPI_Datatype
 * arguments for MPI function calls. The static function, element_type(),
 * returns an argument of type MPI_Datatype that matches a C++ datatype with
 * an MPI_Datatype.
 *
 */
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

template<class T>
struct MPI_Traits 
{
};

//---------------------------------------------------------------------------//
// SPECIALIZATIONS OF MPI_Traits FOR DIFFERENT TYPES
//---------------------------------------------------------------------------//

template<>
struct MPI_Traits<char>
{
    static MPI_Datatype element_type() { return MPI_CHAR; }
};

template<>
struct MPI_Traits<unsigned char>
{
    static MPI_Datatype element_type() { return MPI_UNSIGNED_CHAR; }
};

template<>
struct MPI_Traits<short>
{
    static MPI_Datatype element_type() { return MPI_SHORT; }
};

template<>
struct MPI_Traits<unsigned short>
{
    static MPI_Datatype element_type() { return MPI_UNSIGNED_SHORT; }
};

template<>
struct MPI_Traits<int>
{
    static MPI_Datatype element_type() { return MPI_INT; }
};

template<>
struct MPI_Traits<unsigned int>
{
    static MPI_Datatype element_type() { return MPI_UNSIGNED; }
};

template<>
struct MPI_Traits<long>
{
    static MPI_Datatype element_type() { return MPI_LONG; }
};

template<>
struct MPI_Traits<unsigned long>
{
    static MPI_Datatype element_type() { return MPI_UNSIGNED_LONG; }
};

template<>
struct MPI_Traits<float>
{
    static MPI_Datatype element_type() { return MPI_FLOAT; }
};

template<>
struct MPI_Traits<double>
{
    static MPI_Datatype element_type() { return MPI_DOUBLE; }
};

template<>
struct MPI_Traits<long double>
{
    static MPI_Datatype element_type() { return MPI_LONG_DOUBLE; }
};

} // end namespace rtt_c4

#endif                          // __c4_MPI_Traits_hh__

//---------------------------------------------------------------------------//
//                              end of c4/MPI_Traits.hh
//---------------------------------------------------------------------------//
