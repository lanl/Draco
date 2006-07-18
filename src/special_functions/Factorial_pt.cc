//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   sf/factorial_pt.cc
 * \author Kelly Thompson
 * \date   Mon Nov 8 11:17:12 2004
 * \brief  Provide explicit instantiations of templatized factorial function. 
 * \note   Copyright 2006 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "Factorial.i.hh"

namespace rtt_sf
{

//---------------------------------------------------------------------------//
// Make factorial valid only for int and unsigned.

template
unsigned factorial( unsigned const k ) ;

template
int factorial( int const k ) ;

template
long factorial( long const k ) ;

template
double factorial_fraction( unsigned const k, unsigned const l );

template
double factorial_fraction( int const k, int const l );

template
double factorial_fraction( long const k, long const l );

} // end namespace rtt_sf

//---------------------------------------------------------------------------//
//              end of sf/factorial_pt.cc
//---------------------------------------------------------------------------//
