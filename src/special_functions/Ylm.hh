//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   special_functions/Ylm.hh
 * \author Kent Budge
 * \date   Tue Jul  6 10:03:25 MDT 2004
 * \brief  Declare the Ylm function template.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef special_functions_Ylm_hh
#define special_functions_Ylm_hh

#include "ds++/config.h"

namespace rtt_sf {

/*!
 * \brief Compute the spherical harmonic coefficient multiplied by the appropriate Associated
 *        Legendre Polynomial \f$ c_{l,k}P_{l,k}(\mu) \f$.
 */
double cPlk(unsigned const l, unsigned const m, double const mu);

/*!
 * \brief Compute Morel's Galerkin-quadrature spherical harmonic coefficient multiplied by the
 *        appropriate Associated Legendre Polynomial.
 */
double cPlkGalerkin(unsigned const l, unsigned const m, double const mu, double const sumwt);

//! Compute the normalized spherical harmonic \f$ y_{l,k}(\theta,\phi) \f$.
double normalizedYlk(unsigned const l, int const m, double const theta, double const phi);

//! Compute the real portion of the spherical harmonic \f$ Y_{l,k}(\theta,\phi) \f$.
double realYlk(unsigned const l, int const m, double const theta, double const phi);

//! Compute the imaginary portion of the spherical harmonic \f$ Y_{l,k}(\theta,\phi) \f$.
double complexYlk(unsigned const l, int const m, double const theta, double const phi);

//! Compute the spherical harmonic as used by Morel's Galerkin Quadrature paper.
double galerkinYlk(unsigned const l, int const m, double const mu, double const phi,
                   double const sumwt);

double Ylm(unsigned const l, int const m, double const mu, double const phi, double const sumwt);

} // end namespace rtt_sf

#endif // special_functions_Ylm_hh

//------------------------------------------------------------------------------------------------//
// end of utils/Ylm.hh
//------------------------------------------------------------------------------------------------//
