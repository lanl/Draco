//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   min/linmin.hh
 * \author Kent Budge
 * \date   Tue Aug 17 15:30:23 2004
 * \brief  Find minimum of a multivariate function on a specified line.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

#ifndef min_linmin_hh
#define min_linmin_hh

#include "brent.hh"
#include "mnbrak.hh"
#include <vector>

namespace rtt_min {

//------------------------------------------------------------------------------------------------//
/*!
 * \class f1dim
 * \brief Helper template for line minimization
 *
 * \arg \a RandomContainer A random-access container on double.
 * \arg \a Function A function type supporting <code>double
 *         operator()(RandomContainer const &)</code>.
 */
template <typename RandomContainer, typename Function> class f1dim {
public:
  f1dim(Function func_in, RandomContainer &p_in, RandomContainer &xi_in)
      : func(func_in), p(p_in), xi(xi_in) {}

  double operator()(double const x) const {
    Check(p.size() < UINT_MAX);
    auto n = static_cast<unsigned>(p.size());
    std::vector<double> xt(n);
    for (unsigned i = 0; i < n; ++i) {
      xt[i] = p[i] + x * xi[i];
    }
    return func(xt);
  }

private:
  Function func;
  RandomContainer &p, &xi;
};

//------------------------------------------------------------------------------------------------//
/*!
 * Find minimum of a multivariate function on a specified line.
 *
 * \arg \a RandomContainer A random-access container on double.
 * \arg \a Function A function type supporting <code>double
 *         operator()(RandomContainer const &)</code>.
 *
 * \param[in,out] p  On entry, contains a starting guess of the minimum.  On
 *                   exit, contains an improved estimate of the minimum.
 * \param[in,out] xi On entry, contains the search direction.  On exit,
 *                   contains the actual displacement to the minimum in the
 *                   search direction.
 * \param[out] fret  Final minimum value of the function.
 * \param[in]  func  Function to be minimized
 */

template <typename RandomContainer, typename Function>
void linmin(RandomContainer &p, RandomContainer &xi, double &fret, Function func) {
  using std::numeric_limits;

  double const TOL = sqrt(numeric_limits<double>::epsilon());

  Check(p.size() < UINT_MAX);
  auto const n = static_cast<unsigned>(p.size());

  // Initial guess for brackets
  double ax = 0.0;
  double xx = 1.0;

  double bx, fa, fx, fb, xmin;
  f1dim<RandomContainer, Function> f1(func, p, xi);
  mnbrak<f1dim<RandomContainer, Function>>(ax, xx, bx, fa, fx, fb, f1);

  fret = brent(ax, xx, bx, f1, TOL, xmin);
  for (unsigned j = 0; j < n; ++j) {
    xi[j] *= xmin;
    p[j] += xi[j];
  }
}

} // end namespace rtt_min

#endif // min_linmin_hh

//------------------------------------------------------------------------------------------------//
// end of min/linmin.hh
//------------------------------------------------------------------------------------------------//
