//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   min/powell.hh
 * \author Kent Budge
 * \date   Tue Aug 17 15:30:23 2004
 * \brief  Find minimum of a multivariate function.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

#ifndef min_powell_hh
#define min_powell_hh

#include "linmin.hh"

namespace rtt_min {

//------------------------------------------------------------------------------------------------//
/*!
 * Find minimum of a multivariate function.
 *
 * \arg \a RandomContainer A random-access container on double.
 * \arg \a Function A function type supporting <code>double
 * operator()(RandomContainer const &)</code>.
 *
 * \param[in,out] p On entry, contains a starting guess of the minimum.  On
 * exit, contains an improved estimate of the minimum.
 * \param[in,out] xi On entry, contains the starting set of search directions,
 * which should usually be chosen to be the unit vectors.  On exit, contains
 * the final set of search directions.
 * \param[in] ftol Desired function tolerance.  When a search step fails to
 * decrease the function value by more than ftol, the search is done.
 * \param[in,out] iter Number of iterations to take/taken to complete the
 * search.
 * \param[out] fret Final minimum value of the function.
 * \param[in] func Function to be minimized
 */

template <class RandomContainer, class Function>
void powell(RandomContainer &p, RandomContainer &xi, double const ftol, unsigned &iter,
            double &fret, Function func) {
  using rtt_dsxx::square;
  using std::vector;

  unsigned const ITMAX = iter;
  double const TINY = 1.0e-25;

  Check(p.size() < UINT_MAX);
  auto const n = static_cast<unsigned>(p.size());
  vector<double> pt(n), ptt(n), xit(n);

  fret = func(p);

  std::copy(p.begin(), p.end(), pt.begin());

  double fptt;
  for (iter = 0;; ++iter) {
    double const fp = fret;
    unsigned ibig = 0;
    double del = 0.0;
    for (unsigned i = 0; i < n; ++i) {
      for (unsigned j = 0; j < n; ++j) {
        xit[j] = xi[j + n * i];
      }
      fptt = fret;
      linmin(p, xit, fret, func);
      if (fptt - fret > del) {
        del = fptt - fret;
        ibig = i + 1;
      }
    }
    if (2 * (fp - fret) <= ftol * (fabs(fp) + fabs(fret)) + TINY) {
      return;
    }
    if (iter == ITMAX) {
      throw std::range_error("powell exceeding maximum iterations");
    }
    for (unsigned j = 0; j < n; ++j) {
      ptt[j] = 2 * p[j] - pt[j];
      xit[j] = p[j] - pt[j];
      pt[j] = p[j];
    }
    fptt = func(ptt);
    if (fptt < fp) {
      double const t =
          2 * (fp - 2 * fret + fptt) * square(fp - fret - del) - del * square(fp - fptt);
      if (t < 0.0) {
        linmin(p, xit, fret, func);
        for (unsigned j = 0; j < n; ++j) {
          xi[j + n * (ibig - 1)] = xi[j + n * (n - 1)];
          xi[j + n * (n - 1)] = xit[j];
        }
      }
    }
  }
}

} // end namespace rtt_min

#endif // min_powell_hh

//------------------------------------------------------------------------------------------------//
// end of min/powell.hh
//------------------------------------------------------------------------------------------------//
