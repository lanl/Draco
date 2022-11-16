//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   ds++/DracoMath.hh
 * \author Kent G. Budge
 * \date   Wed Jan 22 15:18:23 MST 2003
 * \brief  New or overloaded cmath or cmath-like functions.
 * \note   Copyright (C) 2013-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_dsxx_DracoMath_hh
#define rtt_dsxx_DracoMath_hh

#include "Constexpr_Functions.hh"
#include "Soft_Equivalence.hh"
#include <algorithm>
#include <complex>
#include <cstdlib>
#include <functional>

namespace rtt_dsxx {

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Return the conjugate of a quantity.
 *
 * The default implementation assumes a field type that is self-conjugate, such as \c double.  An
 * example of a field type that is \em not self-conjugate is \c complex.
 *
 * \tparam Field type
 * \param arg Field type
 */
template <typename Field> inline Field conj(const Field &arg) { return arg; }

// Specializations for non-self-conjugate types
template <> inline std::complex<double> conj(const std::complex<double> &arg) {
  return std::conj(arg);
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Return the cube of a value.
 *
 * \param[in] x Value to be cubed.
 * \return x^3
 *
 * \c Semigroup is a type representing an algebraic structure closed under multiplication such as
 * the integers or the reals.
 */
template <typename Semigroup> inline Semigroup cube(Semigroup const &x) { return x * x * x; }

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Return the positive difference of the arguments.
 *
 * This is a replacement for the FORTRAN DIM function.
 *
 * \arg \a Ordered_Group_Element A type for which operator< and unary operator- are defined and
 *      which can be constructed from a literal \c 0.
 *
 * \param a Minuend
 * \param b Subtrahend
 * \return max(0, a-b)
 *
 * \deprecated A FORTRAN relic that should disappear eventually.
 */
template <typename Ordered_Group_Element>
/* [[deprecated]] */ inline Ordered_Group_Element dim(Ordered_Group_Element a,
                                                      Ordered_Group_Element b) {
  if (a < b)
    return Ordered_Group_Element(0);
  else
    return a - b;
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Return the square of a value.
 *
 * \arg \a Semigroup A type representing an algebraic structure closed under multiplication, such as
 *      the integers or the reals.
 *
 * \param[in] x Value to be squared.
 * \return x^2
 */
template <typename Semigroup> inline Semigroup square(const Semigroup &x) { return x * x; }

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Compute the hypotenuse of a right triangle.
 *
 * This function evaluates the expression \f$\sqrt{a^2+b^2}\f$ in a way that is insensitive to
 * roundoff and preserves range.
 *
 * \arg \a Real A real number type
 * \param a First leg of triangle
 * \param b Second leg of triangle
 * \return Hypotenuse, \f$\sqrt{a^2+b^2}\f$
 */
template <typename Real> inline double pythag(Real a, Real b) {
  using std::abs;
  Real absa = abs(a), absb = abs(b);
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
 * \arg \a Ordered_Group A type for which \c operator< and unary \c operator- are defined and which
 * can be compared to literal \c 0.
 *
 * \param a Argument supplying magnitude of result.
 * \param b Argument supplying sign of result.
 * \return Value |a| with the sign of b.
 */
template <typename Ordered_Group> inline Ordered_Group sign(Ordered_Group a, Ordered_Group b) {
  using std::abs; // just to be clear

  if (b < 0)
    return -abs(a);
  else
    return abs(a);
}
//------------------------------------------------------------------------------------------------//
/*!
 * \brief Fast ceiling of an integer division
 *
 * \param[in] n numerator
 * \param[in] d denominator
 * \return ceil(n/d)
 *
 * If concerned about overflow, consider 1 + ((n-1)/d).
 */
template <typename T>
constexpr inline typename std::enable_if<std::is_integral<T>::value, T>::type
ceil_int_division(T const n, T const d) {
  Require(d != 0);
  return (n + d - 1) / d;
}

} // namespace rtt_dsxx

#endif // rtt_dsxx_DracoMath_hh

//------------------------------------------------------------------------------------------------//
// end of DracoMath.hh
//------------------------------------------------------------------------------------------------//
