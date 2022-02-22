//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   special_functions/F_eta.hh
 * \author Kent Budge
 * \date   Mon Sep 20 15:01:53 2004
 * \brief  For a fermionic species, calculate the dimensionless number density from the
 *         dimensionless chemical potential and dimensionless temperature.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved.
 *
 * These routines are based on C routines from \em Numerical \em Recipes.
 */
//------------------------------------------------------------------------------------------------//

#ifndef sf_F_eta_hh
#define sf_F_eta_hh

#include "ds++/config.h"

namespace rtt_sf {

//! Calculate the relativistic Fermi-Dirac dimensionless number density.
double F_eta(double const eta, double const gamma);

} // end namespace rtt_sf

#endif // sf_F_eta_hh

//------------------------------------------------------------------------------------------------//
// end of sf/F_eta.hh
//------------------------------------------------------------------------------------------------//
