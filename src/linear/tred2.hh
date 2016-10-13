//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   linear/tred2.hh
 * \author Kent Budge
 * \date   Thu Sep  2 14:49:55 2004
 * \brief  Householder reduce a symmetric matrix.
 * \note   Copyright (C) 2016 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef linear_tred2_hh
#define linear_tred2_hh

namespace rtt_linear {

// Use explicit instantiation model.

//! Householder-reduce a symmetric matrix
template <class FieldVector1, class FieldVector2, class FieldVector3>
void tred2(FieldVector1 &a, unsigned n, FieldVector2 &d, FieldVector3 &e);

} // end namespace rtt_linear

#endif // linear_tred2_hh

//---------------------------------------------------------------------------//
// end of linear/tred2.hh
//---------------------------------------------------------------------------//
