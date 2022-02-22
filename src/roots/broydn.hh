//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   roots/broydn.hh
 * \author Kent Budge
 * \date   Wed Jul  7 09:14:09 2004
 * \brief  Find a solution of a set of nonlinear equations.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef roots_broydn_hh
#define roots_broydn_hh

#include "fdjac.hh"
#include "lnsrch.hh"
#include "linear/qr_unpack.hh"
#include "linear/qrdcmp.hh"
#include "linear/qrupdt.hh"
#include "linear/rsolv.hh"
#include "linear/svbksb.hh"
#include "linear/svdcmp.hh"
#include <numeric>

namespace rtt_roots {
//------------------------------------------------------------------------------------------------//
/*!
 * \brief Use Broyden's method to solve a set of nonlinear equations.
 *
 * This procedure searches for a root of a system of nonlinear equations using line minimization
 * along a set of conjugate gradient directions.  When a line minimization fails to make significant
 * progress, the conjugate gradient calculation is restarted.
 *
 * The Jacobian matrix is initially inverted using QR decomposition.  Subsequent Jacobian inverses
 * are estimated from the results of previous minimizations.  If this process breaks down, the
 * procedure goes back to QR decomposition.  If the Jacobian is singular, the singular value
 * decomposition is computed and used to try to get the algorithm out of a tight spot.  This makes
 * for a fairly robust and efficient method.
 *
 * \arg \a Field A field type, such as double or complex.
 * \arg \a Function_N_to_N A multifunctor type, representing a set of N functions in N variables.
 *         Each function returns the residual of the corresponding equation in the nonlinear system.
 *         \sa Function_N_to_N
 *
 * \param x Initial estimate of the solution of the set of equations.  On return, contains the best
 *         solution found.
 * \param vecfunc A Function_N_to_N object representing the set of nonlinear equations.
 * \param alf Success determination parameter for line search.  A value of 0 means that any
 *         reduction in the function value is considered a successful search.
 *
 * \pre \c x.size()>0
 */
/* (no doxygen)
 * \param STPMX Set size parameter.  A large value dials a large initial step in line minimization;
 *         a small value dials a small initial step.  Larger is better unless this takes the
 *         argument to the function outside the function's domain.  A typical choice for this
 *         parameter is 100.
 *
 * \bug KGB: STPMX is not a very useful tuning parameter.  In general, the search parameters are not
 *      very well thought out for this procedure.
 */
template <class Field, class Function_N_to_N>
void broydn(std::vector<Field> &x, const double /*STPMX*/, const Function_N_to_N &vecfunc,
            const double alf) {
  Require(x.size() > 0);

  using std::numeric_limits;
  using std::range_error;
  using std::vector;
  using namespace rtt_linear;
  using namespace rtt_roots;

  Check(x.size() < UINT_MAX);
  const auto n = static_cast<unsigned>(x.size());

  vector<Field> c(n);
  vector<Field> d(n);
  vector<Field> fvcold(n);
  vector<Field> g(n);
  vector<Field> p(n);
  vector<Field> qt(n * n);
  vector<Field> r(n * n);
  vector<Field> rr(n * n);
  vector<Field> s(n);
  vector<Field> t(n);
  vector<Field> w(n);
  vector<Field> v(n * n);
  vector<Field> xold(n);
  vector<Field> fvec(n);

  Field f = fnorm(x, fvec, vecfunc);
  bool restrt = true;
  bool singular = false;
  for (;;) {
    if (restrt || singular) {
      fdjac(x, fvec, r, vecfunc);
      for (unsigned i = 0; i < n; ++i) {
        for (unsigned j = 0; j < n; ++j) {
          rr[i + n * j] = r[i + n * j];
        }
      }
      if (qrdcmp(r, n, c, d)) {
        singular = true;
        svdcmp(rr, n, n, w, v);
      } else {
        singular = false;
        qr_unpack(r, n, c, d, qt);
      }
    } else {
      // Calculate change in x over last step.
      for (unsigned i = 0; i < n; i++)
        s[i] = x[i] - xold[i];
      // Calculate secant.
      for (unsigned i = 0; i < n; i++) {
        double sum = r[i + n * i] * s[i];
        for (unsigned j = i + 1; j < n; j++)
          sum += r[i + n * j] * s[j];
        t[i] = sum;
      }
      bool noisy = true;
      for (unsigned i = 0; i < n; i++) {
        double sum = qt[0 + n * i] * t[0];
        for (unsigned j = 1; j < n; j++)
          sum += qt[j + n * i] * t[j];
        w[i] = fvec[i] - fvcold[i] - sum;
        // As we approach the root, the change in f will begin to be swamped by roundoff noise.
        // Filter out all w that are likely to be noisy.  If all w are noisy, don't try to update
        // the Jacobian.
        if (std::abs(w[i]) >
            numeric_limits<double>::epsilon() * (std::abs(fvec[i]) + std::abs(fvcold[i]))) {
          noisy = false; // this w is not yet swamped by roundoff
        } else {
          w[i] = 0.0; // this w is swamped with noise; leave it out
        }
      }
      if (!noisy) {
        for (unsigned i = 0; i < n; i++) {
          double sum = qt[i + n * 0] * w[0];
          for (unsigned j = 1; j < n; j++)
            sum += qt[i + n * j] * w[j];
          t[i] = sum;
        }
        double scale = std::abs(s[0]); // To avoid overflow
        for (unsigned i = 1; i < n; i++) {
          double const fs = std::abs(s[i]);
          if (fs > scale) {
            scale = fs;
          }
        }
        Check(scale > std::numeric_limits<double>::min());
        // Shouldn't happen, as a negligible change in x should already have triggered a successful
        // return.
        double const rscale = 1 / scale;
        double sum = rtt_dsxx::square(s[0] * rscale);
        for (unsigned i = 1; i < n; i++)
          sum += rtt_dsxx::square(s[i] * rscale);
        double const rnorm2 = 1 / ((scale * sum) * scale);
        // The ordering of the above expression is important to avoid overflow.
        for (unsigned i = 0; i < n; i++)
          s[i] *= rnorm2;
        qrupdt(r, qt, n, t, s);
        // Check singularity.
        for (unsigned i = 0; i < n; i++) {
          if (std::abs(r[i + n * i]) < std::numeric_limits<double>::min())
            throw range_error("broydn: singular Jacobian matrix (1)");
        }
      }
    }
    double fold = f;
    if (!singular) {
      for (unsigned i = 0; i < n; i++) {
        double sum = qt[i + n * 0] * fvec[0];
        for (unsigned j = 1; j < n; j++)
          sum += qt[i + n * j] * fvec[j];
        g[i] = sum;
      }
      for (unsigned i = n - 1; i < n; i--)
      // Should be safe: Require(x.size()>0) --> n>0
      {
        double sum = r[0 + n * i] * g[0];
        for (unsigned j = 1; j <= i; j++)
          sum += r[j + n * i] * g[j];
        g[i] = sum;
      }
      for (unsigned i = 0; i < n; i++) {
        xold[i] = x[i];
        fvcold[i] = fvec[i];
      }
      for (unsigned i = 0; i < n; i++) {
        double sum = qt[i + n * 0] * fvec[0];
        for (unsigned j = 1; j < n; j++)
          sum += qt[i + n * j] * fvec[j];
        p[i] = -sum;
      }
      rsolv(r, n, p);
    } else {
      double wmax = std::abs(w[0]);
      for (unsigned i = 1; i < n; i++) {
        if (w[i] > wmax) {
          wmax = w[i];
        }
      }
      for (unsigned i = 0; i < n; i++) {
        if (std::abs(w[i]) < wmax * numeric_limits<double>::epsilon()) {
          w[i] = 0.0;
        }
        xold[i] = x[i];
        fvcold[i] = fvec[i];
        g[i] = fvec[i];
        p[i] = -fvec[i];
      }
      svbksb(rr, w, v, n, n, p, p);
      restrt = true;
    }
    bool check;
    lnsrch(xold, fold, g, p, x, f, check, fvec, vecfunc, alf, 0.0);
    if (check) {
      for (unsigned i = 0; i < n; i++)
        x[i] = xold[i];
      f = fold;
    }
    if (check) {
      if (restrt) {
        return; // limit of accuracy reached
      } else {
        double test = 0.0;
        double den;
        if (0.5 * n > f) {
          den = 0.5 * n;
        } else {
          den = f;
        }
        for (unsigned i = 0; i < n; i++) {
          double fx = std::abs(x[i]);
          double ff = (fx > 1.0 ? fx : 1.0);
          double const temp = std::abs(g[i]) * ff / den;
          if (temp > test)
            test = temp;
        }
        restrt = true;
      }
    } else {
      restrt = false;
    }
  }
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Use Broyden's method to solve a set of nonlinear equations.
 *
 * This procedure searches for a root of a system of nonlinear equations using line minimization
 * along a set of conjugate gradient directions.  When a line minimization fails to make significant
 * progress, the conjugate gradient calculation is restarted.
 *
 * The Jacobian matrix is initially inverted using QR decomposition.  Subsequent Jacobian inverses
 * are estimated from the results of previous minimizations.  If this process breaks down, the
 * procedure goes back to QR decomposition.  If the Jacobian is singular, the singular value
 * decomposition is computed and used to try to get the algorithm out of a tight spot.  This makes
 * for a fairly robust and efficient method.
 *
 * This variant of the Broyden method assumes that the Jacobian is available analytically.
 *
 * \arg \a Field A field type, such as double or complex.
 *
 * \arg \a Function_N_to_N A multifunctor type, representing a set of N functions in N variables.
 *         Each function returns the residual of the corresponding equation in the nonlinear system.
 *         \sa Function_N_to_N
 *
 * \arg \a Function_N_to_NN A multifunctor type, representing the Jacobian of a set of N functions
 *         in N variables.  Each function returns the residual of the corresponding equation in the
 *         nonlinear system.  \sa Function_N_to_NN
 *
 * \param x Initial estimate of the solution of the set of equations.  On return, contains the best
 *         solution found.

 * \param vecfunc A Function_N_to_N object representing the set of nonlinear equations.
 *
 * \param dvecfunc A Function_N_to_NN object representing the Jacobian of the set of nonlinear
 *         equations.
 *
 * \param alf Success determination parameter for line search.  A value of 0 means that any
 *         reduction in the function value is considered a successful search.
 *
 * \param min_lambda Mimimum line search parameter at which to give up.
 */
/* (no doxygen)
 * \param STPMX Set size parameter.  A large value dials a large initial step in line minimization;
 *         a small value dials a small initial step.  Larger is better unless this takes the
 *         argument to the function outside the function's domain.  A typical choice for this
 *         parameter is 100.
 *
 * \bug KGB: STPMX is not a very useful tuning parameter.  In general, the search parameters are not
 *      very well thought out for this procedure.
 */
template <class Field, class Function_N_to_N, class Function_N_to_NN>
void broydn(std::vector<Field> &x, const double /*STPMX*/, Function_N_to_N vecfunc,
            Function_N_to_NN dvecfunc, const double alf, double const min_lambda) {
  using std::numeric_limits;
  using std::range_error;
  using std::vector;
  using namespace rtt_linear;
  using namespace rtt_roots;

  Check(x.size() < UINT_MAX);
  const auto n = static_cast<unsigned>(x.size());

  vector<Field> c(n);
  vector<Field> d(n);
  vector<Field> fvcold(n);
  vector<Field> g(n);
  vector<Field> p(n);
  vector<Field> qt(n * n);
  vector<Field> r(n * n);
  vector<Field> rr(n * n);
  vector<Field> s(n);
  vector<Field> t(n);
  vector<Field> w(n);
  vector<Field> v(n * n);
  vector<Field> xold(n);
  vector<Field> fvec(n);

  Field f = fnorm(x, fvec, vecfunc);
  bool restrt = true;
  bool singular = false;
  for (;;) {
    if (restrt || singular) {
      dvecfunc(x, fvec, r);
      for (unsigned i = 0; i < n; ++i) {
        for (unsigned j = 0; j < n; ++j) {
          rr[i + n * j] = r[i + n * j];
        }
      }
      if (qrdcmp(r, n, c, d)) {
        singular = true;
        svdcmp(rr, n, n, w, v);
      } else {
        singular = false;
        qr_unpack(r, n, c, d, qt);
      }
    } else {
      // Calculate change in x over last step.
      for (unsigned i = 0; i < n; i++)
        s[i] = x[i] - xold[i];
      // Calculate secant.
      for (unsigned i = 0; i < n; i++) {
        double sum = r[i + n * i] * s[i];
        for (unsigned j = i + 1; j < n; j++)
          sum += r[i + n * j] * s[j];
        t[i] = sum;
      }
      bool noisy = true;
      for (unsigned i = 0; i < n; i++) {
        double sum = qt[0 + n * i] * t[0];
        for (unsigned j = 1; j < n; j++)
          sum += qt[j + n * i] * t[j];
        w[i] = fvec[i] - fvcold[i] - sum;
        // As we approach the root, the change in f will begin to be swamped by roundoff noise.
        // Filter out all w that are likely to be noisy.  If all w are noisy, don't try to update
        // the Jacobian.
        if (std::abs(w[i]) >
            numeric_limits<double>::epsilon() * (std::abs(fvec[i]) + std::abs(fvcold[i]))) {
          noisy = false; // this w is not yet swamped by roundoff
        } else {
          w[i] = 0.0; // this w is swamped with noise; leave it out
        }
      }
      if (!noisy) {
        for (unsigned i = 0; i < n; i++) {
          double sum = qt[i + n * 0] * w[0];
          for (unsigned j = 1; j < n; j++)
            sum += qt[i + n * j] * w[j];
          t[i] = sum;
        }
        double scale = std::abs(s[0]); // To avoid overflow
        for (unsigned i = 1; i < n; i++)
          scale = std::max(scale, std::abs(s[i]));
        Check(scale > std::numeric_limits<double>::min());
        // Shouldn't happen, as a negligible change in x should already have triggered a successful
        // return.
        double const rscale = 1 / scale;
        double sum = rtt_dsxx::square(s[0] * rscale);
        for (unsigned i = 1; i < n; i++)
          sum += rtt_dsxx::square(s[i] * rscale);
        double const rnorm2 = 1 / ((scale * sum) * scale);
        // The ordering of the above expression is important to avoid overflow.
        for (unsigned i = 0; i < n; i++)
          s[i] *= rnorm2;
        qrupdt(r, qt, n, t, s);
        // Check singularity.
        for (unsigned i = 0; i < n; i++) {
          if (std::abs(r[i + n * i]) < std::numeric_limits<double>::min())
            throw range_error("broydn: singular Jacobian matrix (2)");
        }
      }
    }
    double fold = f;
    if (!singular) {
      for (unsigned i = 0; i < n; i++) {
        double sum = qt[i + n * 0] * fvec[0];
        for (unsigned j = 1; j < n; j++)
          sum += qt[i + n * j] * fvec[j];
        g[i] = sum;
      }
      for (int i = n - 1; i >= 0; i--) {
        double sum = r[0 + n * i] * g[0];
        for (int j = 1; j <= i; j++)
          sum += r[j + n * i] * g[j];
        g[i] = sum;
      }
      for (unsigned i = 0; i < n; i++) {
        xold[i] = x[i];
        fvcold[i] = fvec[i];
      }
      for (unsigned i = 0; i < n; i++) {
        double sum = qt[i + n * 0] * fvec[0];
        for (unsigned j = 1; j < n; j++)
          sum += qt[i + n * j] * fvec[j];
        p[i] = -sum;
      }
      rsolv(r, n, p);
    } else {
      double wmax = std::abs(w[0]);
      for (unsigned i = 1; i < n; i++) {
        wmax = std::max(wmax, w[i]);
      }
      for (unsigned i = 0; i < n; i++) {
        if (std::abs(w[i]) < wmax * numeric_limits<double>::epsilon()) {
          w[i] = 0.0;
        }
        xold[i] = x[i];
        fvcold[i] = fvec[i];
        g[i] = fvec[i];
        p[i] = -fvec[i];
      }
      svbksb(rr, w, v, n, n, p, p);
      restrt = true;
    }
    bool check;
    lnsrch(xold, fold, g, p, x, f, check, fvec, vecfunc, alf, min_lambda);
    if (check) {
      for (unsigned i = 0; i < n; i++)
        x[i] = xold[i];
      f = fold;
    }
    if (check) {
      if (restrt) {
        return; // limit of accuracy reached
      } else {
        double test = 0.0;
        double den = std::max(f, 0.5 * n);
        for (unsigned i = 0; i < n; i++) {
          double const temp = std::abs(g[i]) * std::max(std::abs(x[i]), 1.0) / den;
          if (temp > test)
            test = temp;
        }
        restrt = true;
      }
    } else {
      restrt = false;
    }
  }
}

} // end namespace rtt_roots

#endif // roots_broydn_hh

//------------------------------------------------------------------------------------------------//
// end of roots/broydn.hh
//------------------------------------------------------------------------------------------------//
