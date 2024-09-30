//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   special_functions/ExpInt.cc
 * \author Paul Talbot
 * \date   Tue July 26 13:39:00 2011
 * \brief  Iterative methods to calculate Ei(x), E_n(x)
 * \note   Copyright (C) 2011-2024 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "ExpInt.hh"
#include "dsxx/Assert.hh"
#include "dsxx/Soft_Equivalence.hh"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>

namespace rtt_sf {

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Compute general exponential integrals of type Ei(x) or E_n(x).
 *
 * \param[in] n Used to specify order of E_n(x).
 * \param[in] x Used to specify argument of E_n(x).
 * \return E_n(x) evaluated at the argument x.
 *
 * This routine makes use of those described in Numerical Recipes for C++, 3rd Edition (pp 266-270).
 *
 * E_n(x) is calculated either by special case definitions for n=0 or x=0 with n=0 or 1, or the
 * Lentz algorithm if x>1.0, or the digamma series representation for greater x.
 */
double En(unsigned const n, double const x) {
  using std::numeric_limits;

  Insist(!(x < 0.0 || (std::abs(x) < std::numeric_limits<double>::min() && (n == 0 || n == 1))),
         "bad arguments in En");

  const int MAXIT = 100;
  const double EULER = 0.577215664901533;
  const double EPS = numeric_limits<double>::epsilon();
  const double BIG = numeric_limits<double>::max() * EPS;

  double a, b, c, d, del, fact, h, psi, ans(0.0);

  if (n == 0) {
    return exp(-x) / x; //special case
  }

  int const nm1 = n - 1; // n will always be >= 1

  if (rtt_dsxx::soft_equiv(x, 0.0)) {
    return 1.0 / nm1; //also special case
  }

  if (x > 1.0) {
    b = x + n;
    c = BIG;
    d = 1.0 / b;
    h = d;
    for (int i = 1; i <= MAXIT; ++i) {
      a = -1.0 * (i * (nm1 + i));
      b += 2.0;
      d = 1.0 / (a * d + b); //fixes zero denominator
      c = b + a / c;
      del = c * d;
      h *= del;
      if (fabs(del - 1.0) <= EPS) {
        ans = h * exp(-x);
        return ans;
      }
    }
    Insist(false, "continued fraction failed in En"); // LCOV_EXCL_LINE
  } else                                              //do series
  {
    //first term
    if (nm1 != 0)
      ans = 1.0 / nm1;
    else
      ans = -log(x) - EULER;

    fact = 1.0;
    for (int i = 1; i <= MAXIT; ++i) {
      fact *= -x / i;
      if (i != nm1) {
        del = -fact / (i - nm1);
      } else //compute psi(n)
      {
        psi = -EULER;
        for (int ii = 1; ii <= nm1; ii++)
          psi += 1.0 / ii;
        del = fact * (-log(x) + psi);
      }
      ans += del;
      if (fabs(del) < fabs(ans) * EPS)
        return ans;
    }
    Insist(false, "series failed in E_n"); // LCOV_EXCL_LINE
  }

  Insist(false, "Should never get here.");
  return ans;
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Compute general exponential integrals of type Ei(x).
 *
 * \param x Used to specify argument of Ei(x).
 * \return Ei evaluated at the argument x.
 *
 * This routine makes use of those described in Numerical Recipes for C++, 3rd Edition (pp 266-270).
 *
 * Ei is calculated using power series expansion for x<abs(ln(EPS)), where EPS is the required
 * relative error (set to machine precision), as this is the lower limit of the asymptotic series,
 * which is used for greater values of x.  The asymptotic series includes all converging terms until
 * values are less than EPS.
 *
 * Ei(-x) evaluates to the extension -E_1(x).
 */
double Ei(double const x) {
  // Check(x) // x can be any real value.

  using std::numeric_limits;
  const size_t MAXIT = 100;
  const double EULER = 0.57721564901533;
  const double EPS = numeric_limits<double>::epsilon();
  const double FPMIN = numeric_limits<double>::min() / EPS;

  double fact, prev, sum, term;

  if (x <= 0.0) {
    return -1 * En(1, fabs(x));
  }

  if (x < FPMIN)
    return log(x) + EULER; //Avoid failure from underflow

  if (x <= -log(EPS)) //use power series
  {
    sum = 0.0;
    fact = 1.0;
    size_t k(1);
    for (; k <= MAXIT; ++k) {
      fact *= x / static_cast<double>(k);
      term = fact / static_cast<double>(k);
      sum += term;
      if (term < EPS * sum)
        break;
    }
    Insist(k <= MAXIT, "Series failed in Ei");
    return sum + log(x) + EULER;
  } else //Use asymptotic expansion
  {
    sum = 0.0;
    term = 1.0; //starts with second term
    for (size_t k = 1; k <= MAXIT; ++k) {
      prev = term;
      term *= static_cast<double>(k) / x;
      if (term < EPS)
        break; //term approximates relative error

      if (term < prev) {
        sum += term; //still converging
      } else {
        sum -= prev;
        break; //diverging, take off last term
      }
    }
    return exp(x) * (1.0 + sum) / x;
  }
}

} //end namespace rtt_sf

//------------------------------------------------------------------------------------------------//
// end of ExpInt.cc
//------------------------------------------------------------------------------------------------//
