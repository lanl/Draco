//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   roots/lnsrch.i.hh
 * \author Kent Budge
 * \date   Tue Aug 10 13:21:58 2004
 * \brief  Reduce norm of a set of functions on a ray.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef roots_lnsrch_i_hh
#define roots_lnsrch_i_hh

#include "linear/fnorm.hh"

namespace rtt_roots {

//------------------------------------------------------------------------------------------------//
/*! Reduce the norm of a set of functions on a ray.
 *
 * Given a system of equations, a starting point, and a search vector, this procedure reduces of the
 * norm of the set of functions along the ray defined by the starting point and search direction.
 * It is used mainly as an auxiliary to various quasi-Newton methods for solving sets of nonlinear
 * equations.
 *
 * Because this function is primary a helper function for quasi-Newton methods, no great accuracy
 * requirement is placed on the convergence to the minimum.  Since the Newton step is always tried
 * first, we should get rapid convergence just by accepting this step if we are close enough to the
 * root.  If not, it is sufficient that we have *some* improvement, and \c ALF characterizes how
 * much is required.
 *
 * \arg \a RandomContainer A random access container
 * \arg \a Funtion_N_to_N A type representing a set of N functions of N variables.
 *
 * \param xold Starting point for search.
 * \param fold Starting value of the half norm of the set of functions.
 * \param g    Gradient of the set of functions at the starting point.
 * \param p Search direction.  Normally p+xold will be the estimated location of the minimum and
 *          will always be tried first.
 * \param x On exit, contains the new estimate of the minimum location.
 * \param f On exit, contains the new estimate of the minimum norm.
 * \param check On exit, set to true if the search failed.
 * \param fvec Values of the functions at the new minimum point.
 * \param vecfunc Functor whose norm is to be mininized.
 * \param ALF Minimum relative decrease in norm (from starting value) that is acceptable. A
 *          negative value means to accept any search that does not throw an exception.
 * \param min_lambda Value of lambda (line search parameter) at which to give up
 */
template <class RandomContainer, class Function_N_to_N>
void lnsrch(RandomContainer const &xold, double const fold, RandomContainer const &g,
            RandomContainer &p, RandomContainer &x, double &f, bool &check, RandomContainer &fvec,
            Function_N_to_N const &vecfunc, double const ALF, double const min_lambda) {
  Require(g.size() == xold.size());
  Require(p.size() == xold.size());

  using namespace std;
  using rtt_dsxx::square;
  using rtt_linear::fnorm;

  double const eps = std::numeric_limits<typename RandomContainer::value_type>::epsilon();
  const size_t n = xold.size();

  x.resize(n);
  fvec.resize(n);

  // Calculate initial rate of decrease along search direction.  This is used to determine whether a
  // solution is acceptable.
  double slope = g[0] * p[0];
  for (unsigned i = 1; i < n; i++)
    slope += g[i] * p[i];
  if (slope >= 0) {
    check = true;
    return;
  }
  double lambda = 1.0; // Try the Newton step first.
  check = false;
  double lambda_1(0.0), lambda_2(0.0), f1(0.0), f2(0.0);
  while (true) {
    bool lambda_too_small = true;
    for (unsigned i = 0; i < n; i++) {
      x[i] = xold[i] + lambda * p[i];
      if (!rtt_dsxx::soft_equiv(x[i], xold[i], eps))
        lambda_too_small = false;
    }
    if (lambda_too_small || lambda < min_lambda) {
      check = true;
      return;
    }
    try {
      f = fnorm(x, fvec, vecfunc);
      if (ALF < 0.0 || (f < fold + ALF * lambda * slope) || (fold > 0 && f <= 0.5 * fold)) {
        return;
        // Good enough
      } else if (rtt_dsxx::soft_equiv(f, fold, eps)) {
        // Not good enough, and we've stagnated; give up
        check = true;
        return;
      } else {
        // Not good enough; if first try, use quadratic
        if (rtt_dsxx::soft_equiv(lambda, 1.0, eps)) {
          lambda_2 = 0.0;
          lambda_1 = 1.0;
          lambda = max(0.1, -0.5 * slope / (f - fold - slope));
          f2 = fold;
          f1 = f;
        } else {
          // Second and subsequent tries use a cubic.
          lambda_2 = lambda_1;
          lambda_1 = lambda;
          f2 = f1;
          f1 = f;
          double x1 = f1 - slope * lambda_1 - fold;
          double x2 = f2 - slope * lambda_2 - fold;
          double rden = 1 / (lambda_1 - lambda_2);
          double a = rden * (x1 / square(lambda_1) - x2 / square(lambda_2));
          double b = rden * (-lambda_2 * x1 / square(lambda_1) + lambda_1 * x2 / square(lambda_2));
          double det = b * b - 3 * a * slope;
          if (det < 0)
            lambda = 0.5 * lambda;
          else
            lambda = min(0.5 * lambda, max(0.1 * lambda, (sqrt(det) - b) / (3 * a)));
        }
      }
    } catch (...) {
      lambda *= 0.5;
    }
  }
}

} // end namespace rtt_roots

#endif // roots_lnsrch_i_hh

//------------------------------------------------------------------------------------------------//
// end of roots/lnsrch.i.hh
//------------------------------------------------------------------------------------------------//
