//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   min/mrqmin.t.hh
 * \author Kent Budge
 * \date   Fri Aug 7 11:11:31 MDT 2009
 * \brief  Implementation of mrqmin
 * \note   Copyright (C) 2009-2021 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef min_mrqmin_t_hh
#define min_mrqmin_t_hh

#include "ds++/Assert.hh"
#include "ds++/DracoMath.hh"
#include "ds++/dbc.hh"
#include "linear/gaussj.hh"
#include "min/mrqmin.hh"
#include <cmath>

namespace rtt_min {

using namespace std;

//------------------------------------------------------------------------------------------------//
//! Helper function for mrqmin
template <class RandomContainer, class RandomBoolContainer>
void covsrt(RandomContainer &covar, RandomBoolContainer &ia, unsigned const ma,
            unsigned const mfit) {
  for (unsigned i = mfit; i < ma; ++i) {
    for (unsigned j = 0; j < i + 1; ++j) {
      covar[i + ma * j] = covar[j + ma * i] = 0;
    }
  }
  int k = mfit - 1;
  for (int j = ma - 1; j >= 0; --j) {
    if (ia[j]) {
      for (unsigned i = 0; i < ma; ++i) {
        swap(covar[i + ma * k], covar[i + ma * j]);
      }
      for (unsigned i = 0; i < ma; ++i) {
        swap(covar[k + ma * i], covar[j + ma * i]);
      }
      --k;
    }
  }
}

//------------------------------------------------------------------------------------------------//
//! Helper function for mrqmin

template <class RandomContainer, class RandomBoolContainer, class ModelFunction>
void mrqcof(RandomContainer const &x, RandomContainer const &y, RandomContainer const &sig,
            unsigned const ndata, unsigned const m, RandomContainer &a, RandomBoolContainer &ia,
            RandomContainer &alpha, RandomContainer &beta, unsigned const ma, double &chisq,
            ModelFunction funcs) {
  vector<double> dyda(ma);
  unsigned mfit = 0;
  for (unsigned j = 0; j < ma; ++j) {
    if (ia[j]) {
      mfit++;
    }
  }
  for (unsigned j = 0; j < mfit; ++j) {
    for (unsigned k = 0; k <= j; ++k) {
      alpha[j + ma * k] = 0;
      beta[j] = 0;
    }
  }
  chisq = 0.0;
  vector<double> xx(m);
  for (unsigned i = 0; i < ndata; ++i) {
    for (unsigned j = 0; j < m; ++j) {
      xx[j] = x[j + m * i];
    }
    double ymod;
    funcs(xx, a, ymod, dyda);
    double sig2i = 1.0 / (sig[i] * sig[i]);
    double dy = y[i] - ymod;
    unsigned j, l;
    for (j = 0, l = 0; l < ma; ++l) {
      if (ia[l]) {
        double wt = dyda[l] * sig2i;
        unsigned k;
        unsigned mm;
        for (k = 0, mm = 0; mm < l + 1; ++mm) {
          if (ia[mm]) {
            alpha[j + ma * k++] += wt * dyda[mm];
          }
        }
        beta[j++] += dy * wt;
      }
    }
    chisq += dy * dy * sig2i;
  }
  for (unsigned j = 1; j < mfit; ++j) {
    for (unsigned k = 0; k < j; ++k) {
      alpha[k + ma * j] = alpha[j + ma * k];
    }
  }
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Perform a nonlinear least squares fit of data to a model function.
 *
 * \tparam RandomContainer A random access container
 * \tparam RandomBoolContainer A random access bool container
 * \tparam ModelFunction A model function with effective signature:
 *
 * <code>
 *            void funcs(const double,
 *                       vector<double> const &,
 *                       double &,
 *                       vector<double> &)
 * </code>
 *
 * \param x Ordinates of the data points.
 * \param y Values of the data points.
 * \param sig Uncertainty of the data points.
 * \param n Number of data points
 * \param m Number of ordinates
 * \param a Model parameter values
 * \param ia Map of parameters that are actually allowed to vary
 * \param covar Covariance matrix of the final fit
 * \param alpha Curvature matrix
 * \param ma Number of parameters
 * \param chisq Final chi square of the fit
 * \param funcs Model function
 * \param alamda Fit parameter. If less than zero, initialize the algorithm. If
 *             zero, finalize the algorith. Any other value is assumed to be the
 *             value that was returned by a previous iteration.
 */
template <class RandomContainer, class RandomBoolContainer, class ModelFunction>
void mrqmin(RandomContainer const &x, RandomContainer const &y, RandomContainer const &sig,
            unsigned const n, unsigned const m, RandomContainer &a, RandomBoolContainer &ia,
            RandomContainer &covar, RandomContainer &alpha, unsigned const ma, double &chisq,
            ModelFunction funcs, double &alamda) {
  static unsigned mfit;
  static vector<double> oneda;
  static vector<double> atry, beta, da;
  static double ochisq;

  using rtt_linear::gaussj;

  if (alamda < 0.0) {
    atry.resize(ma);
    beta.resize(ma);
    da.resize(ma);
    covar.resize(ma * ma);
    alpha.resize(ma * ma);
    mfit = 0;
    for (unsigned j = 0; j < ma; ++j) {
      if (ia[j]) {
        ++mfit;
      }
    }
    oneda.resize(mfit);
    alamda = 0.001;
    mrqcof(x, y, sig, n, m, a, ia, alpha, beta, ma, chisq, funcs);
    ochisq = chisq;
    for (unsigned j = 0; j < ma; ++j) {
      atry[j] = a[j];
    }
  }
  vector<double> temp(mfit * mfit);
  for (unsigned j = 0; j < mfit; ++j) {
    for (unsigned k = 0; k < mfit; ++k) {
      covar[j + ma * k] = alpha[j + ma * k];
    }
    covar[j + ma * j] = alpha[j + ma * j] * (1 + alamda);
    for (unsigned k = 0; k < mfit; ++k) {
      temp[j + mfit * k] = covar[j + ma * k];
      oneda[j] = beta[j];
    }
  }
  gaussj(temp, mfit, oneda, 1U);
  for (unsigned j = 0; j < mfit; ++j) {
    for (unsigned k = 0; k < mfit; ++k) {
      covar[j + ma * k] = temp[j + mfit * k];
    }
    da[j] = oneda[j];
  }
  if (std::abs(alamda) < std::numeric_limits<double>::min()) {
    covsrt(covar, ia, ma, mfit);
    covsrt(alpha, ia, ma, mfit);
    return;
  }
  unsigned j;
  unsigned l;
  for (j = 0, l = 0; l < ma; ++l) {
    if (ia[l]) {
      atry[l] = a[l] + da[j++];
    }
  }
  mrqcof(x, y, sig, n, m, atry, ia, covar, da, ma, chisq, funcs);
  if (chisq < ochisq) {
    alamda *= 0.1;
    ochisq = chisq;
    for (unsigned jj = 0; jj < mfit; ++jj) {
      for (unsigned k = 0; k < mfit; ++k) {
        alpha[jj + ma * k] = covar[jj + ma * k];
      }
      beta[jj] = da[jj];
    }
    for (unsigned ll = 0; ll < ma; ++ll) {
      a[ll] = atry[ll];
    }
  } else {
    alamda *= 10.0;
    chisq = ochisq;
  }
}

} // end namespace rtt_min

#endif // min_mrqmin_t_hh

//------------------------------------------------------------------------------------------------//
// end of min/mrqmin.t.hh
//------------------------------------------------------------------------------------------------//
