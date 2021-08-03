//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   kde/kde.i.hh
 * \author Mathew Cleveland
 * \date   Nov. 10th 2020
 * \brief  Member definitions of class kde
 * \note   Copyright (C) 2018-2020 Triad National Security, LLC.
 *         All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef kde_kde_i_hh
#define kde_kde_i_hh

#include "ds++/dbc.hh"

namespace rtt_kde {

//------------------------------------------------------------------------------------------------//
/*!
 * \brief
 * epan_kernel basis function used during reconstruction
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
 * \brief
 * log_transform pushes the distribution to log space which is less sensitive to sharp derivatives.
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
 * \brief
 * log_inv_transform bring the distribution back to the original space.
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

//! Lambda to calculate a vector
auto calc_vec = [](const auto &v1, const auto &v2) {
  return std::array<double, 3>{v2[0] - v1[0], v2[1] - v1[1], v2[2] - v1[2]};
};

//! Lambda to calculate vector magnitude
auto calc_mag = [](const auto &v) { return sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]); };

//! Lambda to calculate unit vector
auto calc_unit_vec = [](const auto &v) {
  const double mag = calc_mag(v);
  return std::array<double, 3>{v[0] / mag, v[1] / mag, v[2] / mag};
};

//------------------------------------------------------------------------------------------------//
/*!
 * \brief
 * Calculate the radius given a sphere center and the current location.
 *
 *
 * \param[in] sphere_center the center location (x,y,z) or (r,z) of the sphere
 * \param[in] location data location (x,y,z) or (r,z)
 *
 * \return radius from cell center
 *
 * Test of kde.
 */

inline double kde::calc_radius(const std::array<double, 3> &sphere_center,
                               const std::array<double, 3> &location) const {
  return calc_mag(calc_vec(sphere_center, location));
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief
 * Calculate the arch length between two points (infinitely extended from sphere center) at a specified radius.
 *
 *
 * \param[in] sphere_center the center location (x,y,z) or (r,z) of the sphere
 * \param[in] radius from sphere center to calculate the arch length
 * \param[in] location_1 data location (x,y,z) or (r,z)
 * \param[in] location_2 data location (x,y,z) or (r,z)
 *
 * \return arch length
 *
 * Test of kde.
 */

inline double kde::calc_arch_length(const std::array<double, 3> &sphere_center, const double radius,
                                    const std::array<double, 3> &location_1,
                                    const std::array<double, 3> &location_2) const {
  const std::array<double, 3> v1{calc_unit_vec(calc_vec(sphere_center, location_1))};
  const std::array<double, 3> v2{calc_unit_vec(calc_vec(sphere_center, location_2))};
  const double cos_theta = v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
  return radius * acos(cos_theta);
}
} // end namespace  rtt_kde

#endif // kde_kde_i_hh

//------------------------------------------------------------------------------------------------//
// end of <pkg>/kde.i.hh
//------------------------------------------------------------------------------------------------//
