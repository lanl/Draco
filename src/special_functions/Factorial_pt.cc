//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   special_functions/Factorial_pt.cc
 * \author Kelly Thompson
 * \date   Mon Nov 8 11:17:12 2004
 * \brief  Provide explicit instantiations of templatized factorial function.
 * \note   Copyright (C) 2016-2020 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "Factorial.t.hh"

namespace rtt_sf {

//------------------------------------------------------------------------------------------------//
// Make factorial valid only for int and unsigned.

template unsigned factorial(unsigned const k);
template int factorial(int const k);
template long factorial(long const k);
template double factorial_fraction(unsigned const k, unsigned const l);
template double factorial_fraction(int const k, int const l);
template double factorial_fraction(long const k, long const l);

} // end namespace rtt_sf

//------------------------------------------------------------------------------------------------//
// end of sf/factorial_pt.cc
//------------------------------------------------------------------------------------------------//
