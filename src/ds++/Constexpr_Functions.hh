//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   ds++/Constexpr_Functions.hh
 * \author Alex R. Long
 * \date   Tue July 2 16:56:16 2019
 * \brief  Constexpr versions of math calls
 * \note   Copyright (C) 2019-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//
#ifndef rtt_dsxx_Constexpr_Functions_hh
#define rtt_dsxx_Constexpr_Functions_hh

#include <limits>

namespace rtt_dsxx {

//------------------------------------------------------------------------------------------------//
double constexpr ce_fabs(double const x) {
  if (x >= 0 && x < std::numeric_limits<double>::max())
    return x;
  else if (x < 0)
    return -1.0 * x;
  return std::numeric_limits<double>::quiet_NaN();
}

//------------------------------------------------------------------------------------------------//
//! Simple recursive Newtwon Raphson used in square root
double constexpr sqrtNewtonRaphson(double const x, double const curr, double const prev) {
  return (ce_fabs((curr - prev) / curr) < 1.0e-15)
             ? curr
             : sqrtNewtonRaphson(x, 0.5 * (curr + x / curr), curr);
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Constexpr version of the square root
 *
 * For a finite and non-negative value of "x", returns an approximation for the square root of "x",
 * Otherwise, returns NaN.
 *
 * \note \c HUGE_VAL == \c std::numeric_limits<double>::infinity()
 */
double constexpr ce_sqrt(double const x) {
  return x >= 0 && x < std::numeric_limits<double>::max()
             ? sqrtNewtonRaphson(x, x, 0)
             : std::numeric_limits<double>::quiet_NaN();
}

} // namespace rtt_dsxx

#endif // rtt_dsxx_Constexpr_Functions_hh

//------------------------------------------------------------------------------------------------//
// end of ds++/Constexpr_Functions.hh
//------------------------------------------------------------------------------------------------//
