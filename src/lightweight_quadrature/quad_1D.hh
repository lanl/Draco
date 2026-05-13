//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   lightweight_quadrature/quad_1D.hh
 * \brief  one-dimensional quadratures, including both classical and non-classical quadruple range
 * \note   Copyright (C) 2025-2026 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef QUAD_1D_HH
#define QUAD_1D_HH

#define GSL_RANGE_CHECK_OFF

#include "units/PhysicalConstexprs.hh"
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>

namespace rtt_lightweight_quadrature {

constexpr rtt_units::PhysicalConstexprs<rtt_units::CGS> pc;
auto constexpr cPI = pc.pi();

//------------------------------------------------------------------------------------------------//

// Gauss-Chebyshev
inline auto chebyshev(int N) {
  // calculate quadrature set
  std::vector<double> x(N, 0), w(N, cPI / N);
  for (int n = 0; n < N; n++) {
    x[n] = std::cos((2 * n + 1) * cPI / (2 * N));
  };

  // return quadrature set
  std::vector<std::vector<double>> quad{x, w};
  return quad;
}

// Gauss-Legendre
inline auto gauss(int N) {
  auto J = gsl_matrix_alloc(N, N);
  gsl_matrix_set_zero(J);

  for (int j = 0; j < N - 1; ++j) {
    const double i(j + 1);
    auto beta = std::sqrt(i * i / (4 * i * i - 1));
    gsl_matrix_set(J, j, j + 1, beta);
    gsl_matrix_set(J, j + 1, j, beta);
  }

  gsl_vector *eval = gsl_vector_alloc(N);
  gsl_matrix *evec = gsl_matrix_alloc(N, N);

  gsl_eigen_symmv_workspace *work = gsl_eigen_symmv_alloc(N);
  gsl_eigen_symmv(J, eval, evec, work);
  gsl_eigen_symmv_free(work);
  gsl_matrix_free(J);

  gsl_eigen_symmv_sort(eval, evec, GSL_EIGEN_SORT_VAL_DESC);

  std::vector<double> x(N, 0.0);
  std::vector<double> w(N, 0.0);

  for (int j = 0; j < N; ++j) {
    x[j] = gsl_vector_get(eval, j);
    const double w0 = gsl_matrix_get(evec, 0, j);
    w[j] = 2 * w0 * w0;
  }

  gsl_vector_free(eval);
  gsl_matrix_free(evec);

  std::vector<std::vector<double>> quad{x, w};
  return quad;
}

// double Gauss
inline auto doublegauss(int N) {
  // Gauss-Legendre
  auto leg = gauss(N / 2);
  std::vector<double> const &xx = leg[0];
  std::vector<double> const &wx = leg[1];

  // remap points and weights
  std::vector<double> x(N, 0), w(N, 0);
  for (int n = 0; n < N / 2; n++) {
    x[N / 2 - n - 1] = (1 - xx[n]) / 2;
    w[N / 2 - n - 1] = wx[n] / 2;
  };
  for (int n = N / 2; n < N; n++) {
    x[n] = -x[N - n - 1];
    w[n] = w[N - n - 1];
  };

  // return quadrature set
  std::vector<std::vector<double>> quad{x, w};
  return quad;
}

// see https://gerard-meurant.fr/p_feb29_2008.pdf for Radau and Lobatto with Stieltjes' algorithm

// Gauss-Radau
inline auto radau(int N) {
  auto Jm = gsl_matrix_alloc(N - 1, N - 1);
  gsl_matrix_set_zero(Jm);

  for (int j = 0; j < N - 2; ++j) {
    const double i(j + 1);
    auto beta = std::sqrt(i * i / (4 * i * i - 1));
    gsl_matrix_set(Jm, j, j + 1, beta);
    gsl_matrix_set(Jm, j + 1, j, beta);
    gsl_matrix_set(Jm, j, j, 1.0);
  }
  gsl_matrix_set(Jm, N - 2, N - 2, 1.0);

  // rhs terms
  auto b = gsl_vector_alloc(N - 1);
  gsl_vector_set_zero(b);
  const double betaN = std::sqrt(((N - 1.0) / (2 * N - 3.0)) * ((N - 1.0) / (2 * N - 1.0)));
  gsl_vector_set(b, N - 2, betaN * betaN);

  int s(0);
  gsl_vector *delta = gsl_vector_alloc(N - 1);
  gsl_permutation *p = gsl_permutation_alloc(N - 1);

  gsl_linalg_LU_decomp(Jm, p, &s);
  gsl_linalg_LU_solve(Jm, p, b, delta);

  // construct tridiagonal matrix
  auto J = gsl_matrix_alloc(N, N);
  gsl_matrix_set_zero(J);
  for (int j = 0; j < N - 1; j++) {
    const double i(j + 1);
    auto beta = std::sqrt(i * i / (4 * i * i - 1));
    gsl_matrix_set(J, j, j + 1, beta);
    gsl_matrix_set(J, j + 1, j, beta);
  };

  gsl_matrix_set(J, N - 2, N - 1, betaN);
  gsl_matrix_set(J, N - 1, N - 2, betaN);
  gsl_matrix_set(J, N - 1, N - 1, gsl_vector_get(delta, N - 2) - 1);

  gsl_permutation_free(p);
  gsl_vector_free(delta);

  gsl_vector *eval = gsl_vector_alloc(N);
  gsl_matrix *evec = gsl_matrix_alloc(N, N);

  gsl_eigen_symmv_workspace *work = gsl_eigen_symmv_alloc(N);
  gsl_eigen_symmv(J, eval, evec, work);
  gsl_eigen_symmv_free(work);
  gsl_matrix_free(J);

  gsl_eigen_symmv_sort(eval, evec, GSL_EIGEN_SORT_VAL_DESC);

  std::vector<double> x(N, 0.0);
  std::vector<double> w(N, 0.0);

  for (int j = 0; j < N; ++j) {
    x[j] = gsl_vector_get(eval, j);
    const double w0 = gsl_matrix_get(evec, 0, j);
    w[j] = 2 * w0 * w0;
  }

  gsl_vector_free(eval);
  gsl_matrix_free(evec);

  std::vector<std::vector<double>> quad{x, w};
  return quad;
}

// Double Radau
inline auto doubleradau(int N) {
  // Gauss-Legendre
  auto rad = radau(N / 2);
  std::vector<double> const &x = rad[0];
  std::vector<double> const &w_x = rad[1];

  // remap points and weights
  std::vector<double> mu(N, 0), w(N, 0);
  for (int n = 0; n < N / 2; n++) {
    mu[N / 2 - n - 1] = (1 - x[n]) / 2;
    w[N / 2 - n - 1] = w_x[n] / 2;
  };
  for (int n = N / 2; n < N; n++) {
    mu[n] = -mu[N - n - 1];
    w[n] = w[N - n - 1];
  };

  // return quadrature set
  std::vector<std::vector<double>> quad{mu, w};
  return quad;
}

// Gauss-Lobatto (for completeness only - decisions need to be made to implement in 2D quadrature)
inline auto lobatto(int N) {
  auto Jm_a = gsl_matrix_alloc(N - 1, N - 1);
  auto Jm_b = gsl_matrix_alloc(N - 1, N - 1);
  gsl_matrix_set_zero(Jm_a);
  gsl_matrix_set_zero(Jm_b);

  for (int j = 0; j < N - 2; ++j) {
    const double i(j + 1);
    auto beta = std::sqrt(i * i / (4 * i * i - 1));

    gsl_matrix_set(Jm_a, j, j + 1, beta);
    gsl_matrix_set(Jm_a, j + 1, j, beta);
    gsl_matrix_set(Jm_a, j, j, 1.0);

    gsl_matrix_set(Jm_b, j, j + 1, beta);
    gsl_matrix_set(Jm_b, j + 1, j, beta);
    gsl_matrix_set(Jm_b, j, j, -1.0);
  }
  gsl_matrix_set(Jm_a, N - 2, N - 2, 1.0);
  gsl_matrix_set(Jm_b, N - 2, N - 2, -1.0);

  // rhs terms
  auto b_a = gsl_vector_alloc(N - 1);
  gsl_vector_set_zero(b_a);
  gsl_vector_set(b_a, N - 2, 1.0);

  int s_a(0);
  gsl_vector *delta_a = gsl_vector_alloc(N - 1);
  gsl_permutation *p_a = gsl_permutation_alloc(N - 1);

  gsl_linalg_LU_decomp(Jm_a, p_a, &s_a);
  gsl_linalg_LU_solve(Jm_a, p_a, b_a, delta_a);

  auto b_b = gsl_vector_alloc(N - 1);
  gsl_vector_set_zero(b_b);
  gsl_vector_set(b_b, N - 2, 1.0);

  int s_b(0);
  gsl_vector *delta_b = gsl_vector_alloc(N - 1);
  gsl_permutation *p_b = gsl_permutation_alloc(N - 1);

  gsl_linalg_LU_decomp(Jm_b, p_b, &s_b);
  gsl_linalg_LU_solve(Jm_b, p_b, b_b, delta_b);

  // construct tridiagonal matrix

  auto J = gsl_matrix_alloc(N, N);
  gsl_matrix_set_zero(J);
  for (int j = 0; j < N - 1; j++) {
    const double i(j + 1);
    auto beta = std::sqrt(i * i / (4 * i * i - 1));
    gsl_matrix_set(J, j, j + 1, beta);
    gsl_matrix_set(J, j + 1, j, beta);
  };

  // calculate parameters
  const double alphN = (gsl_vector_get(delta_b, N - 2) + gsl_vector_get(delta_a, N - 2)) /
                       (gsl_vector_get(delta_a, N - 2) - gsl_vector_get(delta_b, N - 2));
  const double betaN =
      std::sqrt(2 / (gsl_vector_get(delta_a, N - 2) - gsl_vector_get(delta_b, N - 2)));

  gsl_matrix_set(J, N - 2, N - 1, betaN);
  gsl_matrix_set(J, N - 1, N - 2, betaN);
  gsl_matrix_set(J, N - 1, N - 1, alphN);

  gsl_permutation_free(p_a);
  gsl_vector_free(delta_a);
  gsl_permutation_free(p_b);
  gsl_vector_free(delta_b);

  gsl_vector *eval = gsl_vector_alloc(N);
  gsl_matrix *evec = gsl_matrix_alloc(N, N);

  gsl_eigen_symmv_workspace *work = gsl_eigen_symmv_alloc(N);
  gsl_eigen_symmv(J, eval, evec, work);
  gsl_eigen_symmv_free(work);
  gsl_matrix_free(J);

  gsl_eigen_symmv_sort(eval, evec, GSL_EIGEN_SORT_VAL_DESC);

  std::vector<double> x(N, 0.0);
  std::vector<double> w(N, 0.0);

  for (int j = 0; j < N; ++j) {
    x[j] = gsl_vector_get(eval, j);
    const double w0 = gsl_matrix_get(evec, 0, j);
    w[j] = 2 * w0 * w0;
  }

  gsl_vector_free(eval);
  gsl_matrix_free(evec);

  std::vector<std::vector<double>> quad{x, w};
  return quad;
}

//------------------------------------------------------------------------------------------------//
// quadratures for orthogonal polynomial generation

/* 
inline auto tanh_sinh(double h, double tol, double a, double b) {
  std::vector<double> x;
  std::vector<double> w;

  x.reserve(12000);
  w.reserve(12000);

  auto pt = [&a, &b](auto &&k, auto &&hd) {
    const double skh = std::tanh(cPI / 2 * std::sinh(k * hd));
    return 0.5 * (b - a) * skh + 0.5 * (b + a);
  };

  auto wt = [&a, &b](auto &&k, auto &&hd) {
    const double dkh = std::cosh(cPI / 2 * std::sinh(k * hd));
    return 0.5 * (b - a) * cPI / 2 * std::cosh(k * hd) / dkh / dkh;
  };

  int k = 0;
  x.push_back(pt(k, h));

  double weight = wt(k, h);
  w.push_back(weight);

  while (weight > tol) {
    ++k;

    // push back point and weight
    x.push_back(pt(k, h));
    x.push_back(pt(-k, h));

    weight = wt(k, h);
    w.push_back(weight);
    w.push_back(weight);
  };

  std::vector<std::vector<double>> quad{x, w};
  return quad;
}
*/

inline auto err_func(double h, double tol, double a, double b) {
  std::vector<double> x;
  std::vector<double> w;

  x.reserve(12000);
  w.reserve(12000);

  auto pt = [](auto &&k, auto &&hd, auto &&ad, auto &&bd) {
    return 0.5 * ((bd - ad) * std::erf(hd * k) + (bd + ad));
  };
  auto wt = [](auto &&k, auto &&hd, auto &&ad, auto &&bd) {
    return (bd - ad) * std::exp(-hd * hd * k * k) / std::sqrt(cPI);
  };

  int k = 0;
  x.push_back(pt(k, h, a, b));

  double weight = wt(k, h, a, b);
  w.push_back(weight);

  while (weight > tol) {
    ++k;

    // push back point and weight
    x.push_back(pt(k, h, a, b));
    x.push_back(pt(-k, h, a, b));

    weight = wt(k, h, a, b);
    w.push_back(weight);
    w.push_back(weight);
  };

  // return quadrature set
  std::vector<std::vector<double>> quad{x, w};
  return quad;
}

//------------------------------------------------------------------------------------------------//
// general quadrature with weight function W
//------------------------------------------------------------------------------------------------//

template <typename F> auto quadruplerange(int N, F W, double a, double b, bool mapped = false) {
  // tanh-sinh quadrature
  const double h = 0.001;
  const double tol = 1.0e-14;
  //auto quad = tanh_sinh(h, tol, a, b);
  auto quad = err_func(h, tol, a, b);

  const size_t M = quad[0].size();

  // coefficient vectors
  std::vector<double> a_k(N + 1, 0);
  std::vector<double> b_k(N + 1, 0);
  std::vector<double> c_k(N + 1, 0);

  b_k[0] = 1;

  // recurrence parameters
  std::vector<double> A(N, 0), B(N - 1, 0);
  for (auto n = 1; n < N + 1; n++) {
    double x_pn_pn = 0;
    double pn_pn = 0;
    double pm_pm = 0;

    // numerical integration
    for (size_t i = 0; i < M; i++) {
      const double x = mapped ? W(quad[0][i]) : quad[0][i];
      auto wt = quad[1][i] * (mapped ? 1 : W(quad[0][i]));

      // calculate legendre polynomials at quadrature points
      double P_n = 0;
      double P_m = 0;
      for (int k = 0; k < n + 1; k++) {
        P_n += b_k[k] * std::pow(x, k);
        P_m += c_k[k] * std::pow(x, k);
      };

      // calculate inner products
      x_pn_pn += P_n * P_n * wt * x;
      pn_pn += P_n * P_n * wt;
      pm_pm += P_m * P_m * wt;
    };

    //printf("iprod(%d) %.12e\n",n,pn_pn);

    if (n == 1) {
      // calculate parameter
      A[0] = x_pn_pn / pn_pn;

      // Gram-Schmidt
      a_k[0] = -A[0] * b_k[0];
      a_k[1] = b_k[0];
      c_k = b_k;
      b_k = a_k;
    } else {
      // calculate parameters
      A[n - 1] = x_pn_pn / pn_pn;
      B[n - 2] = pn_pn / pm_pm;

      // Gram-Schmidt
      a_k[0] = -A[n - 1] * b_k[0] - B[n - 2] * c_k[0];
      for (int k = 1; k < n + 1; k++) {
        a_k[k] = b_k[k - 1] - A[n - 1] * b_k[k] - B[n - 2] * c_k[k];
      };
      c_k = b_k;
      b_k = a_k;
    };
  };

  // Stieltjes (Wilf)
  auto J = gsl_matrix_alloc(N, N);
  gsl_matrix_set_zero(J);

  for (int j = 0; j < N - 1; ++j) {
    auto ax = A[j];
    gsl_matrix_set(J, j, j, ax);

    auto bx = std::sqrt(B[j]);
    gsl_matrix_set(J, j, j + 1, bx);
    gsl_matrix_set(J, j + 1, j, bx);
  }
  auto ax = A[N - 1];
  gsl_matrix_set(J, N - 1, N - 1, ax);

  // solve eigenvalue problem
  gsl_vector *eval = gsl_vector_alloc(N);
  gsl_matrix *evec = gsl_matrix_alloc(N, N);

  gsl_eigen_symmv_workspace *work = gsl_eigen_symmv_alloc(N);
  gsl_eigen_symmv(J, eval, evec, work);
  gsl_eigen_symmv_free(work);
  gsl_eigen_symmv_sort(eval, evec, GSL_EIGEN_SORT_VAL_ASC);

  gsl_matrix_free(J);

  std::vector<double> x(N, 0.0);
  std::vector<double> w(N, 0.0);

  for (int j = 0; j < N; ++j) {
    x[j] = gsl_vector_get(eval, j);
    const double w0 = gsl_matrix_get(evec, 0, j);
    w[j] = w0 * w0;
  }

  // return quadrature set
  std::vector<std::vector<double>> qr{x, w};
  return qr;
}

//------------------------------------------------------------------------------------------------//
// quadruple range
//------------------------------------------------------------------------------------------------//

//------------------------------------------------------------------------------------------------//
// azimuthal quadruple range (Spence, symmetric about 90 degrees)
inline auto az_qrs90(int N) {
  // quadruple range (remapped from [0,1] to [0,2*pi])
  auto wf = [](auto x) { return 0.5 * std::sin(x); };
  auto qr = quadruplerange(N / 2, wf, 0.0, cPI / 2, true);

  // calculate points and weights
  std::vector<double> x(N, 0), w(N, 0);
  for (int n = 0; n < N / 2; n++) {
    const double beta = std::asin(2 * qr[0][n]);
    const double wt = cPI / 2 * qr[1][n];

    x[n] = std::cos(beta);
    w[n] = wt;
    x[N - n - 1] = -x[n];
    w[N - n - 1] = wt;
  };

  // return quadrature set
  std::vector<std::vector<double>> quad{x, w};
  return quad;
}

//------------------------------------------------------------------------------------------------//
// azimuthal quadruple range (Abu-Shumays, symmetric about 90 degrees)
// note that the points are (Pi/2 - beta_i) where beta_i are the az_qrs90 points
inline auto az2D_qra90(int N) {
  // quadruple range (remapped from [0,1] to [0,2*pi])
  auto wf = [](auto x) { return std::sin(x); };
  auto qr = quadruplerange(N / 2, wf, 0.0, cPI / 2, true);

  // calculate points and weights
  std::vector<double> x(N, 0), w(N, 0);

  for (int n = 0; n < N / 2; n++) {
    const double beta = qr[0][n];
    const double wt = cPI / 2 * qr[1][n];

    x[N / 2 - n - 1] = beta;
    w[N / 2 - n - 1] = wt;
    x[n + N / 2] = -beta;
    w[n + N / 2] = wt;
  };

  // return quadrature set
  std::vector<std::vector<double>> quad{x, w};
  return quad;
}

//------------------------------------------------------------------------------------------------//
// azimuthal quadruple range (Spence, symmetric about 45 degrees)
inline auto az_qrs45(int N) {
  // quadruple range
  constexpr double b = 1.0 / std::numbers::sqrt2;
  auto wf = [](auto x) { return 1.0 / std::sqrt(1 - x * x); };
  auto qr = quadruplerange(N / 2, wf, -b, b);

  // calculate points and weights
  std::vector<double> x(N, 0), w(N, 0);
  for (int n = 0; n < N / 2; n++) {
    const double beta = std::asin(qr[0][n]) + cPI / 4;
    const double wt = cPI / 2 * qr[1][n];

    x[n] = std::cos(beta);
    w[n] = wt;

    x[n + N / 2] = std::cos(beta + cPI / 2);
    w[n + N / 2] = wt;
  }

  // return quadrature set
  std::vector<std::vector<double>> quad{x, w};
  return quad;
}

//------------------------------------------------------------------------------------------------//
// azimuthal quadruple range (Abu-Shumays, symmetric about 45 degrees)
inline auto az_qra45(int N) {
  // quadruple range
  const double b = std::sin(cPI / 8);
  auto wf = [](auto x) { return 2.0 / std::sqrt(1 - x * x); };
  auto qr = quadruplerange(N / 2, wf, -b, b);

  // calculate points and weights
  std::vector<double> x(N, 0), w(N, 0);
  for (int n = 0; n < N / 2; n++) {
    const double beta = 2 * std::asin(qr[0][n]) + cPI / 4;
    const double wt = cPI / 2 * qr[1][n];

    x[n] = std::cos(beta);
    w[n] = wt;

    x[n + N / 2] = std::cos(beta + cPI / 2);
    w[n + N / 2] = wt;
  }

  // return quadrature set
  std::vector<std::vector<double>> quad{x, w};
  return quad;
}

//------------------------------------------------------------------------------------------------//
// azimuthal hextuple range (Spence, symmetric about 30 degrees)
inline auto az_qrs30(int N) {
  if (N % 3) {
    printf("quadrature requires N be a multiple of 3 but N is %d\n", N); // LCOV_EXCL_LINE
    throw;                                                               // LCOV_EXCL_LINE
  }

  // quadruple range
  const double b = 0.5;
  auto wf = [](auto x) { return 1.0 / std::sqrt(1 - x * x); };
  auto qr = quadruplerange(N / 3, wf, -b, b);

  // calculate points and weights
  std::vector<double> x(N, 0), w(N, 0);
  for (int n = 0; n < N / 3; n++) {
    const double beta = std::asin(qr[0][n]) + cPI / 6;
    const double wt = cPI / 3 * qr[1][n];

    x[n] = std::cos(beta);
    w[n] = wt;
    x[n + N / 3] = std::cos(beta + cPI / 3);
    w[n + N / 3] = wt;
    x[n + 2 * N / 3] = std::cos(beta + 2 * cPI / 3);
    w[n + 2 * N / 3] = wt;
  };

  // return quadrature set
  std::vector<std::vector<double>> quad{x, w};
  return quad;
}

//------------------------------------------------------------------------------------------------//
// azimuthal hextuple range (Abu-Shumays, symmetric about 30 degrees)
inline auto az_qra30(int N) {
  if (N % 3) {
    printf("quadrature requires N be a multiple of 3 but N is %d\n", N); // LCOV_EXCL_LINE
    throw;                                                               // LCOV_EXCL_LINE
  }

  // quadruple range
  const double b = std::sin(cPI / 12);
  auto wf = [](auto x) { return 1.0 / std::sqrt(1 - x * x); };
  auto qr = quadruplerange(N / 3, wf, -b, b);

  // calculate points and weights
  std::vector<double> x(N, 0), w(N, 0);
  for (int n = 0; n < N / 3; n++) {
    const double beta = 2 * std::asin(qr[0][n]) + cPI / 6;
    const double wt = cPI / 3 * qr[1][n];

    x[n] = std::cos(beta);
    w[n] = wt;
    x[n + N / 3] = std::cos(beta + cPI / 3);
    w[n + N / 3] = wt;
    x[n + 2 * N / 3] = std::cos(beta + 2 * cPI / 3);
    w[n + 2 * N / 3] = wt;
  };

  // return quadrature set
  std::vector<std::vector<double>> quad{x, w};
  return quad;
}

//------------------------------------------------------------------------------------------------//
// polar quadruple range (Spence, symmetric about 45 degrees)
// note that polar-QRA45 is identical to azimuthal-QRA45
inline auto pol_qrs45(int N) {
  // quadruple range
  const double b = 1 / std::numbers::sqrt2;
  auto wf = [](auto x) {
    return (x + std::sqrt(1 - x * x)) / std::numbers::sqrt2 / std::sqrt(1 - x * x);
  };
  auto qr = quadruplerange(N / 2, wf, -b, b);

  // calculate points and weights
  std::vector<double> x(N, 0), w(N, 0);
  for (int n = 0; n < N / 2; n++) {
    const double beta = std::asin(qr[0][n]) + cPI / 4;
    const double wt = qr[1][n];

    x[n] = std::cos(beta);
    w[n] = wt;
    x[N - n - 1] = -x[n];
    w[N - n - 1] = wt;
  };

  // return quadrature set
  std::vector<std::vector<double>> quad{x, w};
  return quad;
}

//------------------------------------------------------------------------------------------------//
// polar quadruple range (Spence, symmetric about 90 degrees)
// note Spence calls this Polar-QRS90b, while Polar-QRS90a is double-Gauss
// also note that this is also Polar2d_QRA90 (or Abu-Shumays' polar set suitable
// for 2D rectangular geometries)
inline auto pol_qrs90(int N) {
  // quadruple range
  auto wf = [](auto x) { return std::sqrt(1 - x * x); };
  auto qr = quadruplerange(N / 2, wf, 0.0, 1.0, true);

  // calculate points and weights
  std::vector<double> x(N, 0), w(N, 0);
  for (int n = 0; n < N / 2; n++) {
    const double beta = std::asin(qr[0][n]);
    const double wt = qr[1][n];

    x[n] = std::cos(beta);
    w[n] = wt;
    x[N - n - 1] = -x[n];
    w[N - n - 1] = wt;
  };

  //return quadrature set
  std::vector<std::vector<double>> quad{x, w};
  return quad;
}

} // end namespace rtt_lightweight_quadrature

#endif // QUAD_1D_HH
