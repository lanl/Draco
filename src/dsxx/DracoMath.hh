//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   dsxx/DracoMath.hh
 * \author Kent G. Budge
 * \date   Wed Jan 22 15:18:23 MST 2003
 * \brief  New or overloaded constexpr-friendly cmath or cmath-like functions.
 * \note   Copyright (C) 2024-2026 Triad National Security, LLC., All rights reserved.
 *
 * A header that adds lightweight, constexpr‑friendly math utilities to Draco: overloads for conj,
 * cube, square (with optional overflow checks), pythag, sign, and ceil_int_division, plus related
 * helper includes. */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_dsxx_DracoMath_hh
#define rtt_dsxx_DracoMath_hh

#include "Constexpr_Functions.hh"
#include "Soft_Equivalence.hh"
#include "dsxx/device_config.h"
#include <algorithm>
#include <complex>
#include <concepts>
#include <cstdlib>
#include <functional>

namespace rtt_dsxx {

//! A concept type is only defined if multiplication of type T is valid.
template <typename T>
concept Multipliable = requires(T x) { x *x; };

/*! Ordered_Group A concept type for which `operator<` and unary `operator-` are defined and which
 *  can be compared to literal `0`. */
template <typename T>
concept OrderedGroup = requires(T a) {
  { a < T{0} } -> std::convertible_to<bool>;
  { -a } -> std::same_as<T>;
  { std::abs(a) } -> std::same_as<T>;
};

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Return the conjugate of a quantity.
 *
 * The default implementation assumes a field type that is self-conjugate, such as `double`.  An
 * example of a field type that is \em not self-conjugate is `complex`.
 *
 * \tparam Field type
 * \param arg Field type
 */
template <typename Field> constexpr Field conj(const Field &arg) { return arg; }

// Specializations for non-self-conjugate types
template <> constexpr std::complex<double> conj(const std::complex<double> &arg) {
  return std::conj(arg);
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Return the cube of a value.
 *
 * \tparam Semigroup is a type representing an algebraic structure closed under multiplication such
 *    as the integers or the reals.
 *
 * \param[in] x Value to be cubed.
 * \return x^3
 *
 * \sa square
 */
template <Multipliable Semigroup>
constexpr auto cube(Semigroup const &x) noexcept(noexcept(x * x * x)) {
  return x * x * x;
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Return the square of a value.
 *
 * \tparam Semigroup A type representing an algebraic structure closed under multiplication, such as
 *    the integers or the reals.
 *
 * \param[in] x Value to be squared.
 * \return x^2
 *
 * Notes:
 * - Faster than pow(x,2) because there is only one instruction, no call overhead, and compile-time
 *   evaluation.
 * - It would be useful to check for overflow, but this would significantly reduce optimization
 *   since the code would no longer be `constexpr`.
 * - Overflow checks are performed when `DRACO_DIAGNOSTICS_LEVEL_2` is defined.
 *
 * \sa cube
 */
template <Multipliable Semigroup>
constexpr auto square(const Semigroup &x) noexcept(noexcept(x * x)) {
#ifdef DRACO_DIAGNOSTICS_LEVEL_2
  if constexpr (std::integral<Semigroup>) {
#if defined(__clang__) || defined(__GNUC__)
    Semigroup r;
    Check(!(__builtin_mul_overflow(x, x, &r))); // Overflow
    return r;
#else
    if constexpr (std::signed_integral<Semigroup>) {
      Check(!(x > 0 && x > std::numeric_limits<Semigroup>::max() / x)); // Integer Overflow
      Check(!(x < 0 && x < std::numeric_limits<Semigroup>::max() / x)); // Integer Overflow
    } else {
      Check(!(x != 0 && x > std::numeric_limits<Semigroup>::max() / x)); // FP Overflow
    }
    return x * x;
#endif
  } else if constexpr (std::floating_point<Semigroup>) {
    Semigroup r = x * x;
    Check(std::isfinite(r));
    return r;
  } else {
    static_assert(std::is_arithmetic_v<Semigroup>,
                  "square<Semigroup>: T must be an integral or floating-point type");
  }
#else
  return x * x;
#endif
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Compute the hypotenuse of a right triangle.
 *
 * This function evaluates the expression \f$\sqrt{a^2+b^2}\f$ in a way that is insensitive to
 * roundoff and preserves range.
 *
 * \tparam Real A real number type
 * \param a First leg of triangle
 * \param b Second leg of triangle
 * \return Hypotenuse, \f$\sqrt{a^2+b^2}\f$
 */
template <std::floating_point Real> constexpr double pythag(Real a, Real b) {
  Real absa = std::abs(a);
  Real absb = std::abs(b);
  // We must avoid (a/b)^2 > max.
  if (absa <= absb * std::sqrt(std::numeric_limits<Real>::min()))
    return absb;
  if (absb <= absa * std::sqrt(std::numeric_limits<Real>::min()))
    return absa;
  // The regular case...
  if (absa > absb)
    return absa * std::sqrt(1.0 + square(absb / absa));
  else
    return absb * std::sqrt(1.0 + square(absa / absb));
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief  Transfer the sign of the second argument to the first argument.
 *
 * This is a replacement for the FORTRAN SIGN function.  It is useful in numerical algorithms that
 * are roundoff or overflow insensitive and should not be deprecated.
 *
 * \param a Argument supplying magnitude of result.
 * \param b Argument supplying sign of result.
 * \return Value |a| with the sign of b.
 */
template <OrderedGroup T> GPU_HOST_DEVICE constexpr T sign(T a, T b) {
  return (b < 0) ? -std::abs(a) : std::abs(a);
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Fast ceiling of an integer division; works for negative and positive integers.
 *
 * \param[in] n numerator
 * \param[in] d denominator
 * \return ceil(n/d)
 */
template <std::integral T> constexpr T ceil_int_division(T n, T d) {
  Require(d != 0);
  T q = n / d;
  T r = n % d;
  return q + (r != 0 && ((r > 0) == (d > 0)));
}

} // namespace rtt_dsxx

#endif // rtt_dsxx_DracoMath_hh

//------------------------------------------------------------------------------------------------//
// end of DracoMath.hh
//------------------------------------------------------------------------------------------------//
