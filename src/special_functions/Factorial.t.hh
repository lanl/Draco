//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   special_functions/Factorial.t.hh
 * \author Kelly Thompson
 * \date   Mon Nov 8 11:17:12 2004
 * \brief  Provide implementation of templatized factorial function.
 * \note   Copyright (C) 2016-2021 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef sf_factorial_i_hh
#define sf_factorial_i_hh

#include "ds++/Assert.hh"
#include "special_functions/Factorial.hh"
#include <array>

namespace rtt_sf {

//------------------------------------------------------------------------------------------------//
/*!
 * \brief factorial
 *
 * This is a recursive function that computes k!
 *
 * \param k A number whose factorial is desired.
 * \return \f$n!\f$
 * \post \c Result>=1
 */
template <typename T> T factorial(T const k) {
  size_t constexpr N(13);
  // only initialize this once (keyword: static)
  std::array<T, N> constexpr tabularValue = {
      1,    1,     2,      6,       24,       120,      720,
      5040, 40320, 362880, 3628800, 39916800, 479001600
      // These are the next two values in the series.  Unfortunately, they are too big to be held by
      // type long.
      // , 6227020800, 87178291200
  };

  if (k <= 1)
    return 1;
  if (static_cast<unsigned int>(k) < N)
    return tabularValue[k];
  else
    throw std::range_error("factorial out of range");
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Fraction of factorials, \f$ (k!)/(l!) \f$
 *
 * \param k The numerator, \f$ k! \f$
 * \param l The denominator, \f$ l! \f$
 * \return \f$ (k!)/(l!) \f$
 * \post \c Result>=1
 */
template <typename T> double factorial_fraction(T const k, T const l) {
  double result(1.0);
  if (k > l)
    for (T i = l + 1; i <= k; ++i)
      result *= static_cast<double>(i);
  if (k < l)
    for (T i = k + 1; i <= l; ++i)
      result /= static_cast<double>(i);
  return result;
}

} // end namespace rtt_sf

#endif // sf_factorial_t_hh

//------------------------------------------------------------------------------------------------//
// end of sf/factorial.t.hh
//------------------------------------------------------------------------------------------------//
