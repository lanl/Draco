//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   roots/zbrac.i.hh
 * \author Kent Budge
 * \date   Tue Aug 17 15:30:23 2004
 * \brief  Bracket a root of a function.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef roots_zbrac_i_hh
#define roots_zbrac_i_hh

#include "ds++/DracoMath.hh"

namespace rtt_roots {

//------------------------------------------------------------------------------------------------//
/*!
 * Bracket a root of a function.
 *
 * On entry, <code>[x1,x2]</code> defines an initial search interval.  On successful return, the
 * search interval <code>[x1,x2]</code> is guaranteed to contain an odd number of roots if \c func
 * is analytic in \c <code>[x1,x2]</code>.
 *
 * \arg \a Function The template type must support <code>double operator()(double)</code>.
 *
 * \param[in] func Function whose root is to be bracketed
 * \param[in,out] x1 Lower end of search interval
 * \param[in,out] x2 Upper end of search interval
 *
 * \pre <code>func!=NULL</code>
 * \pre <code>x2>x1</code>
 * \post <code>func(x1)<=0 && func(x2)>=0 || func(x1)>=0 && func(x2)<=0 </code>
 *
 * \throw std::domain_error Indicates that no root could be bracketed, perhaps because the function
 * has none.
 *
 * \note zbrac cannot fail to bracket a root or roots of an odd polynomial.  It is otherwise
 * impossible to guarantee success. It can only be guaranteed that, if zbrac returns rather than
 * throwing std::domain_error, then the interval <code>[x1,x2]</code> contains at least one root or
 * singularity of the function.
 */
template <typename Function, typename Real> void zbrac(Function func, Real &x1, Real &x2) {
  Require(x2 > x1);

  using namespace std;

  Real const eps = std::numeric_limits<Real>::epsilon();

  Real f1 = func(x1);
  Real f2 = func(x2);
  Real f3;
  Real af1 = (f1 > 0. ? f1 : -f1);
  Real af2 = (f2 > 0. ? f2 : -f2);

  double scale = 0.5;
  while ((f1 < 0 && f2 < 0) || (f1 > 0 && f2 > 0)) {
    if (af1 < af2) {
      Real x0 = x1 - scale * (x2 - x1);
      if (rtt_dsxx::soft_equiv(x0, x1, eps)) {
        throw std::domain_error("zbrac: could not find search interval (1)");
      }
      try {
        Real f0 = func(x0);
        if (!std::isnormal(f0))
          throw std::runtime_error("");
        f1 = f0;
        af1 = (f1 > 0. ? f1 : -f1);
        x1 = x0;
        if (scale < 0.5)
          scale *= 1.5;
      } catch (...) {
        scale *= 0.5;
      }
    } else if (af1 > af2) {
      Real const x3 = x2 + scale * (x2 - x1);
      // Leaving this in place causes cdi_analytic_tstAnalytic_EoS to fail.
      //
      // if (rtt_dsxx::soft_equiv(x2, x3, eps))
      // if( x2 == x3 ) {
      //   throw std::domain_error("zbrac: could not find search interval");
      // }
      try {
        Real const ff3 = func(x3);
        if (std::isinf(ff3))
          throw std::runtime_error("");
        f2 = ff3;
        af2 = (f2 > 0. ? f2 : -f2);
        x2 = x3;
        if (scale < 0.5)
          scale *= 1.5;
      } catch (...) {
        scale *= 0.5;
      }
    } else {
      Real x0 = x1 - 0.5 * scale * (x2 - x1);
      Real x3 = x2 + scale * (x2 - x1);
      if (rtt_dsxx::soft_equiv(x0, x1, eps) || rtt_dsxx::soft_equiv(x2, x3, eps)) {
        throw std::domain_error("zbrac: could not find search interval (2)");
      }
      try {
        Real f0 = func(x0);
        if (!std::isnormal(f0))
          throw std::runtime_error("");
        f3 = func(x3);
        f1 = f0;
        af1 = (f1 > 0. ? f1 : -f1);
        x1 = x0;
        f2 = f3;
        af2 = (f2 > 0. ? f2 : -f2);
        x2 = x3;
        if (scale < 0.5)
          scale *= 1.5;
      } catch (...) {
        scale *= 0.5;
      }
    }
  }

  Ensure((std::isnormal(f1) && std::isnormal(f2)) ? (f1 <= 0 && f2 >= 0) || (f1 >= 0 && f2 <= 0)
                                                  : true);
}

} // end namespace rtt_roots

#endif // roots_zbrac_i_hh

//------------------------------------------------------------------------------------------------//
// end of roots/zbrac.i.hh
//------------------------------------------------------------------------------------------------//
