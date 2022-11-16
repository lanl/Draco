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

#include <cstddef>
#include <vector>

namespace rtt_dsxx {
namespace interpolate {

//! Perform a 1-d linear interpolation
double linear_1d(double const x1, double const x2, double const y1, double const y2,
                 double const x);

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

} // namespace interpolate
} // namespace rtt_dsxx

#endif
