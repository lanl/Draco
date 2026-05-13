//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   norms/L2norm.hh
 * \author Kent Budge
 * \brief  Define template function L2norm
 * \note   Copyright (C) 2010-2026 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef norms_L2norm_hh
#define norms_L2norm_hh

#include "c4/C4_Functions.hh"
#include "norms/norm.hh"
#include <cmath>
#include <concepts>
#include <iostream>
#include <numeric>
#include <ranges>

//! A concept type that is a container of floating point values. eg. `vector<double>` or
//! `list<float>`.
template <typename R>
concept floating_point_random_access_range =
    std::ranges::random_access_range<R> && std::floating_point<std::ranges::range_value_t<R>>;

namespace rtt_norms {

//------------------------------------------------------------------------------------------------//
//! Helper type for L2norm.
template <std::floating_point T> double accumulate_norm_(T const init, T value) {
  return init + norm(value);
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief compute the L2 norm of the container `x`.
 * \param x Container representing a real vector whose norm is desired.
 */
template <floating_point_random_access_range In> double L2norm(In const &x) {
  double value =
      std::accumulate(x.begin(), x.end(), 0.0, accumulate_norm_<std::ranges::range_value_t<In>>);

  rtt_c4::global_sum(value);
  auto x_length(x.size());

  rtt_c4::global_sum(x_length);
  Require(x_length > 0);

  value = std::sqrt(value / static_cast<double>(x_length));

  Ensure(value >= 0.0);
  return value;
}

//------------------------------------------------------------------------------------------------//
/*!
 * This function computes the norm of the difference between two vectors. We have found that this is
 * a surprisingly common operation, and there is advantage to not having to compute the difference
 * vector if all we want is its norm.
 *
 * \tparam In1 An input container type whose elements are real.
 * \tparam In2 An input container type whose elements are real.
 *
 * \param x Container representing a real vector.
 * \param y Container representing a real vector.
 */
template <floating_point_random_access_range In1, floating_point_random_access_range In2>
double L2norm_diff(In1 const &x, In2 const &y) {
  Require(std::ranges::size(x) == std::ranges::size(y));
  double value = 0.0;

  auto xi = std::ranges::begin(x);
  auto yi = std::ranges::begin(y);

  for (; xi != std::ranges::end(x); ++xi, ++yi) {
    value += norm_diff(*xi, *yi);
  }

  rtt_c4::global_sum(value);

  Check(x.size() < UINT_MAX);
  auto x_length(static_cast<unsigned>(x.size()));

  rtt_c4::global_sum(x_length);
  Require(x_length > 0);

  value = sqrt(value / x_length);

  Ensure(value >= 0.0);
  return value;
}

} // end namespace rtt_norms

#endif // norms_L2norm_hh

//------------------------------------------------------------------------------------------------//
// end of norms/L2norm.hh
//------------------------------------------------------------------------------------------------//
