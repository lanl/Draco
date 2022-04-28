//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   roots/zbrent.i.hh
 * \author Kent Budge
 * \date   Tue Aug 17 15:57:06 2004
 * \brief  Find a bracketed root of a function.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef roots_zbrent_i_hh
#define roots_zbrent_i_hh

#include "ds++/DracoMath.hh"
#include <limits>

namespace rtt_roots {

//------------------------------------------------------------------------------------------------//
/*!
 * Pinpoint a bracketed root of a function.
 *
 * \c zbrent returns when either \c tol or \c ftol achieves the requested value or when \c ITMAX
 *           iterations have been attempted.
 *
 * \arg \a Function The template type must support <code>double operator()(double)</code>.
 *
 * \param[in] func Function whose zero is to be found.
 * \param[in] x1 Lower limit of search range
 * \param[in] x2 Upper limit of search range.  \c x1 and \c x2 must bracket a root; that is, the
 *           value of \c func(x1) must differ in sign from the value of \c func(x2).
 *           rtt_roots::zbrac may be helpful for bracketing a root.
 * \param[in] itmax Maximum number of iterations to try.
 * \param[in,out] tol On entry, the desired absolute tolerance in the argument. On return, this is
 *           replaced with the tolerance actually achieved.
 * \param[in,out] ftol On entry, the desired absolute tolerance in the function value.  On return,
 *           this is replaced with the tolerance actually achieved.
 *
 * \pre <code>(fa>=0.0 && fb<=0.0) || (fa<=0.0 && fb>=0.0)</code>
 *
 * \return An estimated zero of the function.
 *
 * \throw std::domain_error If the function does not appear to be analytic in the search interval.
 *
 * \note The tolerances \c tol and \c ftol should always be checked after the function is called, to
 *       verify that an acceptable accuracy has been achieved.  If the function is analytic in
 *       <code>[x1,x2]</code>, then in principle zbrent cannot fail.  In practice, zbrent may fail
 *       to achieve tight tolerances if the function is either very slowly or very rapidly varying
 *       near its zero, due to roundoff.
 */
template <typename Function, typename Real>
Real zbrent(Function func, Real x1, Real x2, unsigned itmax, Real &tol, Real &ftol) {
  using std::min;
  using std::numeric_limits;

  double const eps = std::numeric_limits<Real>::epsilon();
  Real a = x1, b = x2, c = x2;
  Real fa = func(a), fb = func(b);

  Require((fa >= 0.0 && fb <= 0.0) || (fa <= 0.0 && fb >= 0.0));

  Real fc = fb;
  Real d = b - a, e = b - a, xm = 0;
  for (unsigned iter = 0; iter < itmax; iter++) {
    if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
      c = a;
      fc = fa;
      e = d = b - a;
    }
    Real afc = (fc > 0. ? fc : -fc);
    Real afb = (fb > 0. ? fb : -fb);
    if (afc < afb) {
      a = b;
      b = c;
      c = a;
      fa = fb;
      fb = fc;
      fc = fa;
    }

    Real absb = (b > 0. ? b : -b);
    Real const tol1 = 2.0 * numeric_limits<double>::epsilon() * absb + 0.5 * tol;

    xm = 0.5 * (c - b);
    Real axm = (xm > 0. ? xm : -xm);
    afb = (fb > 0. ? fb : -fb);
    if (axm <= tol1 || afb <= ftol) {
      tol = fabs((xm));
      ftol = fabs((fb));
      return b;
    }
    Real ae = (e > 0. ? e : -e);
    Real afa = (fa > 0. ? fa : -fa);
    if (ae >= tol1 && afa > afb) {
      Real const s = fb / fa;
      Real p, q;
      if (rtt_dsxx::soft_equiv(a, c, eps)) {
        p = 2.0 * xm * s;
        q = 1.0 - s;
      } else {
        q = fa / fc;
        Real const r = fb / fc;
        p = s * (2.0 * xm * q * (q - r) - (b - a) * (r - 1.0));
        q = (q - 1.0) * (r - 1.0) * (s - 1.0);
      }
      if (p > 0.0) {
        q = -q;
      }
      p = (p > 0. ? p : -p);
      Real const tlq = tol1 * q;
      Real const atlq = (tlq > 0. ? tlq : -tlq);
      Real const min1 = 3.0 * xm * q - atlq;
      Real eq = e * q;
      Real const min2 = (eq > 0. ? eq : -eq);
      if (2.0 * p < min(min1, min2)) {
        e = d;
        d = p / q;
      } else {
        d = xm;
        e = d;
      }
    } else {
      d = xm;
      e = d;
    }
    a = b;
    fa = fb;
    Real const ad = (d > 0. ? d : -d);
    if (ad > tol1) {
      b += d;
    } else {
      if (xm < 0.0) {
        b -= tol1;
      } else {
        b += tol1;
      }
    }
    try {
      fb = func(b);
      if (std::isinf(fb)) {
        throw std::domain_error("function is not analytic in the search interval");
      }
    } catch (...) {
      throw std::domain_error("function is not analytic in the search interval");
    }
  }
  tol = fabs(xm);
  ftol = fabs(b);
  return b;
}

} // end namespace rtt_roots

#endif // roots_zbrent_i_hh

//------------------------------------------------------------------------------------------------//
// end of roots/zbrent.i.hh
//------------------------------------------------------------------------------------------------//
