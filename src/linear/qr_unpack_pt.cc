//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   linear/qr_unpack.cc
 * \author Kent Budge
 * \date   Wed Aug 11 15:21:38 2004
 * \brief  Specializations of qr_unpack
 * \note   Copyright (C) 2016-2017 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <vector>

#include "qr_unpack.i.hh"

namespace rtt_linear {
using std::vector;

//---------------------------------------------------------------------------//
// RandomContainer = vector<double>
//---------------------------------------------------------------------------//

template DLL_PUBLIC_linear void qr_unpack(vector<double> &r, const unsigned n,
                                          const vector<double> &c,
                                          const vector<double> &d,
                                          vector<double> &qt);

} // end namespace rtt_linear

//---------------------------------------------------------------------------//
// end of qr_unpack.cc
//---------------------------------------------------------------------------//
