//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   fit/svdfit.i.hh
 * \author Kent Budge
 * \date   Mon Aug  9 13:17:31 2004
 * \brief  Calculate the singular value decomposition of a matrix.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef fit_svdfit_i_hh
#define fit_svdfit_i_hh

#include "ds++/DracoMath.hh"
#include "linear/svbksb.hh"
#include "linear/svdcmp.hh"

namespace rtt_fit {

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Compute a generalized least squares fit.
 *
 * Given a set of data, find the best linear combination of arbitrary basis functions to fit the
 * data.
 *
 * \tparam RandomContainer A random access container type
 * \tparam Functor A functor type taking a double value and a reference to a RandomContainer that
 *         stores the values of the basis functions for the double value in the RandomContainer.
 *
 * \param x Ordinates of data set. For multivariable fits, one can let the ordinates be an index
 *            into a table containing the full vector of independent variables for each data point.
 * \param y Abscissae of data set.
 * \param sig Uncertainty in the data. Where unknown or not applicable, one can set all values to 1.
 * \param a On entry, must be sized to the number of basis functions. On exit, contains the
 *            coefficients of the fit.
 * \param u On exit, contains the upper matrix of the singular value decomposition of the fit.
 * \param v On exit, containts the diagonal matrix of the singular value decomposition of the fit,
 *            e.g., the singular values.
 * \param w On exit, containts the lower matrix of the singular value decomposition of the fit.
 * \param chisq On return, contains the chi square of the fit (meaasure of goodness of fit.)
 * \param funcs Functor to calculate the basis functions for a given argument.
 * \param[in] TOL reset denormalized w-values below TOL*max(w) to a hard-zero.
 */
template <typename RandomContainer, typename Functor>
void svdfit(RandomContainer const &x, RandomContainer const &y, RandomContainer const &sig,
            RandomContainer &a, RandomContainer &u, RandomContainer &v, RandomContainer &w,
            double &chisq, Functor &funcs, double TOL) {
  Require(x.size() == y.size());
  Require(x.size() == sig.size());
  Require(a.size() > 0);

  using rtt_dsxx::square;
  using rtt_linear::svbksb;
  using rtt_linear::svdcmp;

  Check(x.size() < UINT_MAX);
  Check(a.size() < UINT_MAX);
  auto const ndata = static_cast<unsigned>(x.size());
  auto const ma = static_cast<unsigned>(a.size());

  std::vector<double> b(ndata);
  std::vector<double> afunc(ma);

  u.resize(ndata * ma);

  for (unsigned i = 0; i < ndata; ++i) {
    funcs(x[i], afunc);
    double const tmp = 1.0 / sig[i];
    for (unsigned j = 0; j < ma; ++j) {
      u[i + ndata * j] = afunc[j] * tmp;
    }
    b[i] = y[i] * tmp;
  }
  svdcmp(u, ndata, ma, w, v);
  double wmax = 0.0;
  for (unsigned j = 0; j < ma; ++j) {
    if (w[j] > wmax) {
      wmax = w[j];
    }
  }
  double const thresh = TOL * wmax;
  for (unsigned j = 0; j < ma; ++j) {
    if (w[j] < thresh) {
      w[j] = 0.0;
    }
  }
  svbksb(u, w, v, ndata, ma, b, a);
  chisq = 0.0;
  for (unsigned i = 0; i < ndata; ++i) {
    funcs(x[i], afunc);
    double sum = 0.0;
    for (unsigned j = 0; j < ma; ++j) {
      sum += a[j] * afunc[j];
    }
    chisq += square((y[i] - sum) / sig[i]);
  }
}

} // end namespace rtt_fit

#endif // fit_svdfit_i_hh

//------------------------------------------------------------------------------------------------//
// end of fit/svdfit.i.hh
//------------------------------------------------------------------------------------------------//
