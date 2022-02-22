//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   special_functions/F2inv.hh
 * \author Kent Budge
 * \date   Tue Sep 21 09:20:10 2004
 * \brief  Inverse Fermi-Dirac integral of 1/2 order.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef sf_F2inv_hh
#define sf_F2inv_hh

#include "ds++/config.h"

namespace rtt_sf {
//! Compute the inverse Fermi-Dirac function of index 2.
double F2inv(double f);

//! Compute the inverse Fermi-Dirac function of index 2 and its derivative.
void F2inv(double f, double &mu, double &dmudf);

} // end namespace rtt_sf

#endif // sf_F2inv_hh

//------------------------------------------------------------------------------------------------//
// end of sf/F2inv.hh
//------------------------------------------------------------------------------------------------//
