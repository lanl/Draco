//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   ds++/DracoMath.hh
 * \author Kent G. Budge
 * \date   Wed Jan 22 15:18:23 MST 2003
 * \brief  New or overloaded cmath or cmath-like functions.
 * \note   Copyright (C) 2016-2021 Triad National Security, LLC., All rights reserved. */
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
// isFinite.hh
//
// Try to use the C++11/C99 functions isinf, isnan and isfinite defined in <cmath> instead of
// defining our own.  I would like to use C++11 implemenations which are true functions in the std::
// namespace.  The problem here is that PGI/13.7 does not have these language features.  However,
// PGI does provide the C99 _macros_ of the same name (w/o namespace qualifier).
//------------------------------------------------------------------------------------------------//
#if defined _WIN32 || defined __CYGWIN__

template <typename T> bool isNan(T a) { return _isnan(a); }
template <typename T> bool isInf(T a) { return !_finite(a); }
template <typename T> bool isFinite(T a) { return _finite(a); }

#elif defined draco_isPGI

template <typename T> bool isNan(T a) { return isnan(a); }
template <typename T> bool isInf(T a) { return isinf(a); }
template <typename T> bool isFinite(T a) { return isfinite(a); }

#else

template <typename T> bool isNan(T a) { return std::isnan(a); }
template <typename T> bool isInf(T a) { return std::isinf(a); }
template <typename T> bool isFinite(T a) { return std::isfinite(a); }

#endif

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
 * \brief Do a linear interpolation between two values.
 *
 * \param[in] x1 x coordinate of first data point.
 * \param[in] y1 y coordinate of first data point.
 * \param[in] x2 x coordinate of second data point.
 * \param[in] y2 y coordinate of second data point.
 * \param[in] x  x coordinate associated with requested y value.
 * \return The y value associated with x based on linear interpolation between (x1,y1) and (x2,y2).
 *
 * Given two points (x1,y1) and (x2,y2), use linaer interpolation to find the y value associated
 * with the provided x value.
 *
 *          y2-y1
 * y = y1 + ----- * (x-x1)
 *          x2-x1
 *
 * \pre  x in (x1,x2), extrapolation is not allowed.
 * \post y in (y1,y2), extrapolation is not allowed.
 */
constexpr inline double linear_interpolate(double const x1, double const x2, double const y1,
                                           double const y2, double const x) {
  Require(ce_fabs(x2 - x1) > std::numeric_limits<double>::epsilon());
  Require(((x >= x1) && (x <= x2)) || ((x >= x2) && (x <= x1)));

  // return value
  double const value = (y2 - y1) / (x2 - x1) * (x - x1) + y1;

  Ensure(((value >= y1) && (value <= y2)) || ((value >= y2) && (value <= y1)));
  return value;
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
