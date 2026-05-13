//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   norms/norm.hh
 * \author Kent Budge
 * \brief  Declare template class norm
 * \note   Copyright (C) 2010-2026 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef norms_norm_hh
#define norms_norm_hh

#include "dsxx/DracoMath.hh"
#include <cmath>
#include <concepts>

namespace rtt_norms {

//------------------------------------------------------------------------------------------------//
//! Compute the norm of a scalar. This is a name alias for rtt_dsxx:square.
constexpr auto norm = [](const auto &x) noexcept(noexcept(rtt_dsxx::square(x))) {
  return rtt_dsxx::square(x);
};

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Compute the norm of the difference of two scalars.
 *
 * \tparam T A floating point type
 *
 * \param[in] x The difference between x and y will be squared.
 * \param[in] y The difference between x and y will be squared.
 * \return `(x-y)^2`
 *
 * \sa square
 */
template <std::floating_point T>
constexpr auto norm_diff(const T &x, const T &y) noexcept(noexcept((x - y) * (x - y))) {
  const T d = x - y;
#ifdef DRACO_DIAGNOSTICS_LEVEL_2
  T r = d * d;
  Check(std::isfinite(r));
  return r;
#else
  return d * d;
#endif
}

} // end namespace rtt_norms

#endif // norms_norm_hh

//------------------------------------------------------------------------------------------------//
// end of norms/norm.hh
//------------------------------------------------------------------------------------------------//
