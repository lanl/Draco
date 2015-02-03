//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   linear/qrdcmp_pt.cc
 * \author Kent Budge
 * \date   Wed Aug 11 15:21:38 2004
 * \brief  Specializations of qrdcmp
 * \note   Copyright (C) 2006-2015 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <vector>
#include "qrdcmp.i.hh"

namespace rtt_linear
{
using std::vector;

//---------------------------------------------------------------------------//
// RandomContainer=vector<double>
//---------------------------------------------------------------------------//

template DLL_PUBLIC
bool qrdcmp(vector<double> &a,
    unsigned n,
    vector<double> &c,
    vector<double> &d);

} // end namespace rtt_linear

//---------------------------------------------------------------------------//
// end of qrdcmp_pt.cc
//---------------------------------------------------------------------------//
