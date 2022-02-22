//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   special_functions/Power.hh
 * \author Mike Buksas
 * \date   Thu Jul 20 17:23:31 2006
 * \brief  A meta-programming implementation of the Russian Pesant algorithm.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved.
 *
 * Use meta-programming to generate an efficient routine to compute integer powers.
 *
 * E.g.: Power<4>(x) computes the fourth power of x. The function constructed at compile time should
 * be equivalent to the following:
 *
 * \code
 *   v1 = x * x
 *  v2 = v1 * v1;
 *  return v2;
 * \endcode
 *
 * Assuming that the compiler is inlining aggressively.
 *
 * Likewise, Power<7>(x) should generate code equivalent to:
 *
 * \code
 *  v1 = x * x
 *  v2 = v1 * x;
 *  v3 = v1 * v2;
 *  return v3;
 * \endcode
 *
 * The meta-algorithm is based on the Russian Peasant Algorithm, modified to be recursive, since
 * this is required for template meta-programming in C++.
 */
//------------------------------------------------------------------------------------------------//

#ifndef special_functions_Power_hh
#define special_functions_Power_hh

namespace rtt_sf {

// Protect the implementation detail of struct P from accidental usage outside of this file.

namespace {

/* Struct P implements a static method: 'compute' which recusvively calls P::compute for N/2.  We
 * use a struct to hold compute because we need to specialize for N=0 and this is not possible with
 * a template function.
 */

template <int N, typename F> struct P {
  static F compute(F x, F p) {
    x *= x;
    constexpr bool is_divby2 = ((N / 2) * 2 == N);
    if (is_divby2)
      return P<N / 2, F>::compute(x, p);
    else
      return P<N / 2, F>::compute(x, x * p);
  }
};

// Specialize struct P on N=0 to terminate the recursion.
template <typename F> struct P<0, F> {
  static F compute(F /*x*/, F p) { return p; }
};
} // namespace

/* Function Power recursively implements the first half of the Russian Pesant algorithm, by
 * repeatedly computing x=x^2, N=N/2, so long as the remaining exponent is even. When an odd
 * exponent is reached, it dispatches to P<N>::compute(x,x) which contiunues the calculation.
 */

template <int N, typename F> F Power(F x) {
  constexpr bool is_zero = (N == 0);
  constexpr bool is_divby2 = ((N / 2) * 2 == N);
  if (is_zero)
    return static_cast<F>(1);
  else if (is_divby2)
    return Power<N / 2>(x * x);
  else
    return P<N / 2, F>::compute(x, x);
}

} // end namespace rtt_sf

#endif // special_functions_Power_hh

//------------------------------------------------------------------------------------------------//
// end of special_functions/Power.hh
//------------------------------------------------------------------------------------------------//
