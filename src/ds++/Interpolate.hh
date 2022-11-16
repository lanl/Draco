//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   ds++/Interpolate.hh
 * \author Kendra P. Long, Kent G. Budge
 * \date   Wed Jan 22 15:18:23 MST 2003
 * \brief  Interpolation functions.
 * \note   Copyright (C) 2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//
#ifndef rtt_compton_tools_Interpolate_hh
#define rtt_compton_tools_Interpolate_hh

#include "ds++/Constexpr_Functions.hh"
#include "ds++/dbc.hh"
#include <cstddef>
#include <limits>
#include <vector>

namespace rtt_dsxx {
namespace interpolate {

//! Perform a 3-d (trilinear) interpolation
double linear_3d(double const x0, double const x1, double const y0, double const y1,
                 double const z0, double const z1, double const f000, double const f100,
                 double const f001, double const f101, double const f010, double const f110,
                 double const f011, double const f111, double const x, double const y,
                 double const z);

//! Compute the Lagrange multipliers for set of interp. regions, each with n_local pointsh
std::vector<double> lagrange_multipliers(const size_t n_break, const size_t n_local,
                                         const std::vector<double> &points);

//! Perform 1-D Lagrange polynomial interpolation, given y values
double lagrange_1d(const std::vector<double> &data, const std::vector<double> &xs,
                   const std::vector<double> &cxs, const double x);

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Do a 1D linear interpolation between two values.
 *
 * \param[in] x1 x coordinate of first data point.
 * \param[in] y1 y coordinate of first data point.
 * \param[in] x2 x coordinate of second data point.
 * \param[in] y2 y coordinate of second data point.
 * \param[in] x  x coordinate associated with requested y value.
 * \return The y value associated with x based on linear interpolation between (x1,y1) and (x2,y2).
 *
 * Given two points (x1,y1) and (x2,y2), use linear interpolation to find the y value associated
 * with the provided x value.
 *
 *          y2-y1
 * y = y1 + ----- * (x-x1)
 *          x2-x1
 *
 * \pre  x in (x1,x2), extrapolation is not allowed.
 * \post y in (y1,y2), extrapolation is not allowed.
 */
constexpr inline double linear_1d(double const x1, double const x2, double const y1,
                                  double const y2, double const x) {
  Require(ce_fabs(x2 - x1) > std::numeric_limits<double>::epsilon());
  Require(((x >= x1) && (x <= x2)) || ((x >= x2) && (x <= x1)));

  // return value
  double const value = (y2 - y1) / (x2 - x1) * (x - x1) + y1;

  Ensure(((value >= y1) && (value <= y2)) || ((value >= y2) && (value <= y1)));
  return value;
}

} // namespace interpolate
} // namespace rtt_dsxx

#endif
