//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   sf/F_eta_inv.hh
 * \author Kent Budge
 * \date   Mon Sep 20 15:01:53 2004
 * \brief  For a fermion species, calculate the dimensionless chemical 
 * potential given the dimensionless number density and dimensionless
 * temperature.  
 * \note   � Copyright 2006 LANSLLC All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef sf_F_eta_inv_hh
#define sf_F_eta_inv_hh

namespace rtt_sf
{

//! Calculate the relativistic Fermi-Dirac dimensionless chemical potential.
double F_eta_inv(double n, double gamma);

} // end namespace rtt_sf

#endif // sf_F_eta_inv_hh

//---------------------------------------------------------------------------//
//              end of sf/F_eta_inv.hh
//---------------------------------------------------------------------------//
