//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   shared_lib/test/Foo.cc
 * \author Rob Lowrie
 * \date   Wed Dec 29 11:24:17 2004
 * \brief  Implementatin of Foo.
 * \note   Copyright (C) 2016 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "Foo.hh"

namespace rtt_shared_lib_test {

Foo::Foo(const double x) : d_base(x) {}

double Foo::compute(const double x) const { return x * x + d_base; }

} // end namespace rtt_shared_lib_test

//---------------------------------------------------------------------------//
// end of shared_lib/test/Foo.cc
//---------------------------------------------------------------------------//
