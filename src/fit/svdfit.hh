//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   fit/svdfit.hh
 * \author Kent Budge
 * \date   Mon Aug  9 13:17:31 2004
 * \brief  Calculate the singular value decomposition of a matrix.
 * \note   Copyright (C) 2010-2026 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef utils_svdfit_hh
#define utils_svdfit_hh

#include "dsxx/DracoMath.hh"
#include "dsxx/dll_declspec.h"
#include "linear/svbksb.hh"
#include "linear/svdcmp.hh"
#include <concepts>
#include <cstddef>
#include <numeric>
#include <ranges>
#include <type_traits>

namespace rtt_fit {

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Floating‑point random‑access container
 *
 *  – The type must model std::random_access_range (so we can use operator[] and know the size).
 *  – The elements stored in the container must be a floating‑point type (float, double, long
 *    double, or any user‑defined type that satisfies std::floating_point). */
//------------------------------------------------------------------------------------------------//
template <typename C>
concept FloatingPointRandomContainer =
    std::ranges::random_access_range<C> &&              // random‑access
    std::floating_point<std::ranges::range_value_t<C>>; // floating-point type

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Compute a generalized least squares fit using singular value decomposition (SVD).
 *
 * Given a set of data, find the best linear combination of arbitrary basis functions to fit the
 * data.
 *
 * \tparam RandomContainer A random access container type, must be of type std::random_access_range
 * \tparam Functor A functor type taking a double value and a reference to a RandomContainer that
 *         stores the values of the basis functions for the double value in the RandomContainer.
 *
 * \param[in] x Ordinates of data set. For multivariable fits, one can let the ordinates be an index
 *            into a table containing the full vector of independent variables for each data point.
 * \param[in] y Abscissae of data set.
 * \param[in] sig Uncertainty in the data. Where unknown or not applicable, one can set all values
 *            to 1.
 * \param a On entry, must be sized to the number of basis functions. On exit, contains the
 *            coefficients of the fit.
 * \param u On exit, contains the upper matrix of the singular value decomposition of the fit.
 * \param v On exit, containts the diagonal matrix of the singular value decomposition of the fit,
 *            e.g., the singular values.
 * \param w On exit, containts the lower matrix of the singular value decomposition of the fit.
 * \param chi_squared On return, contains the chi‑square statistic of the fit (measure of goodness
 *            of fit.)
 * \param funcs Functor to calculate the basis functions for a given argument.
 * \param[in] TOL reset denormalized w-values below TOL*max(w) to a hard-zero.
 *
 * \pre `x.size() == y.size() == sig.size()` and `a.size() > 0`.
 * \pre `sig[i] != 0` for all i (division by zero would be undefined).
 * \pre `funcs` must fill a container of size `a.size()` for each call.
 * \post `a` holds the least‑squares coefficients for the supplied basis functions.
 * \post `u`, `v`, and `w` contain the SVD of the weighted design matrix.
 * \post `chi_squared` equals \f$\sum_i ((y_i - \sum_j a_j f_j(x_i))/sig_i)^2\f$.
 * \post Singular values smaller than `TOL * max(w)` are set to zero, improving numerical stability.
 */
template <FloatingPointRandomContainer RandomContainer, typename Functor>
void svdfit(RandomContainer const &x, RandomContainer const &y, RandomContainer const &sig,
            RandomContainer &a, RandomContainer &u, RandomContainer &v, RandomContainer &w,
            double &chi_squared, Functor const &funcs, double const TOL = 1.0e-13) {
  Require(x.size() == y.size());
  Require(x.size() == sig.size());
  Require(a.size() > 0);
  Require(std::ranges::all_of(sig, [](auto const &s) { // Require sig[i] != 0
    return !rtt_dsxx::soft_equiv(s, 0.0, 1.0e-14);
  }));

  auto const ndata = x.size();
  auto const ma = a.size();

  std::vector<double> b(ndata);
  std::vector<double> afunc(ma);
  u.resize(ndata * ma);

  for (size_t i = 0; i < ndata; ++i) {
    funcs(x[i], afunc);
    double const tmp = 1.0 / sig[i];
    for (size_t j = 0; j < ma; ++j) {
      u[i + ndata * j] = afunc[j] * tmp;
    }
    b[i] = y[i] * tmp;
  }
  rtt_linear::svdcmp(u, static_cast<unsigned>(ndata), static_cast<unsigned>(ma), w, v);
  double const wmax = std::ranges::max(w | std::views::take(ma)); // max element in range.
  double const thresh = TOL * wmax;
  for (double &w_i : w | std::views::take(ma)) { // set epsilon values -> 0.0
    if (w_i < thresh)
      w_i = 0.0;
  }
  rtt_linear::svbksb(u, w, v, static_cast<unsigned>(ndata), static_cast<unsigned>(ma), b, a);
  chi_squared = 0.0;
  for (size_t i = 0; i < ndata; ++i) {
    funcs(x[i], afunc);
    double const sum = std::inner_product(a.begin(), a.begin() + ma, afunc.begin(), 0.0);
    chi_squared += rtt_dsxx::square((y[i] - sum) / sig[i]);
  }
}

} // end namespace rtt_fit

#endif // utils_svdfit_hh

//------------------------------------------------------------------------------------------------//
// end of utils/svdfit.hh
//------------------------------------------------------------------------------------------------//
