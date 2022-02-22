//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   kde/kde.i.hh
 * \author Mathew Cleveland
 * \date   Nov. 10th 2020
 * \brief  Member definitions of class kde
 * \note   Copyright (C) 2021-2022 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

#ifndef kde_kde_i_hh
#define kde_kde_i_hh

#include "ds++/dbc.hh"

namespace rtt_kde {

//------------------------------------------------------------------------------------------------//
/*!
 * \brief epan_kernel 
 *
 * Basis function used during reconstruction.
 *
 * Epanechnikov kenrel to be used in reconstrtuction
 *
 * \param[in] x from kernel origin
 * \return distribution weight based on distance from the kernel center 
 *
 * Test of kde.
 */
inline double kde::epan_kernel(const double x) const {
  const double x2 = x * x;
  return x2 > 1.0 ? 0.0 : 0.75 * (1.0 - x2);
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief log_transform 
 *
 * Transforms data to log space given a bias.
 *
 *
 * \param[in] value of original distribution
 * \param[in] bias used to ensure positivity
 *
 * \return the logarithmic transform of the original value
 *
 * Test of kde.
 */
inline double kde::log_transform(const double value, const double bias) const {
  Require(value + bias > 0.0);
  return log(value + bias);
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief log_inv_transform 
 *
 * Inverse transform back from log space given the current bias.
 *
 *
 * \param[in] log_value of original distribution
 * \param[in] bias used to ensure positivity
 *
 * \return the logarithmic transform of the original value
 *
 * Test of kde.
 */

inline double kde::log_inv_transform(const double log_value, const double bias) const {
  return exp(log_value) - bias;
}

} // end namespace  rtt_kde

#endif // kde_kde_i_hh

//------------------------------------------------------------------------------------------------//
// end of <pkg>/kde.i.hh
//------------------------------------------------------------------------------------------------//
