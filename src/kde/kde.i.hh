//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   kde/kde.i.hh
 * \author Mathew Cleveland
 * \date   Nov. 10th 2020
 * \brief  Member definitions of class kde
 * \note   Copyright (C) 2021-2023 Triad National Security, LLC., All rights reserved.
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
 * \brief int_epan_kernel 
 *
 * Integral form of the basis function used during reconstruction (note this will be zero if the
 * integration area is zero - so don't use in the point wise reconstruction).
 *
 * Epanechnikov kenrel to be used in reconstrtuction
 *
 * \param[in] x0 beginning point of the integration
 * \param[in] x1 ending point of the integration
 * \return distribution weight based on distance from the kernel center 
 *
 * Test of kde.
 */
inline double kde::int_epan_kernel(const double x0, const double x1) const {
  // fix the integration bounds from -1..1
  Require(!(x0 > x1));
  const double a = std::min(std::max(x0, -1.0), 1.0);
  const double b = std::min(std::max(x1, -1.0), 1.0);
  const double a3 = a * a * a;
  const double b3 = b * b * b;
  // int_a^b(3/4(*1.0-x**2)dx = 1/3*3/4*(a**3-3a-b**3+3b)
  const double weight = 0.25 * (a3 - 3.0 * a - b3 + 3.0 * b);
  Ensure(!(weight < 0.0));
  Ensure(!(weight > 1.0));
  return weight;
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
