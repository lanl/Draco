//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   linear/qrdcmp.t.hh
 * \author Kent Budge
 * \date   Mon Aug  9 13:17:31 2004
 * \brief  Calculate the Q-R decomposition of a square matrix.
 * \note   Copyright (C) 2016-2021 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef linear_qrdcmp_t_hh
#define linear_qrdcmp_t_hh

#include "qrdcmp.hh"
#include "ds++/Assert.hh"
#include "ds++/DracoMath.hh"
#include <cmath>

namespace rtt_linear {

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Compute the QR decomposition of a square matrix.
 *
 * Compute the decomposition of a matrix \f$A=QR\f$ where \f$Q\f$ is an orthonormal matrix and
 * \f$R\f$ is upper triangular.  This procedure is specialized for a square matrix and does no
 * pivoting (which is necessary only if the matrix is very close to singular.)  The QR decomposition
 * is about twice as expensive as the LU decomposition and so should be used only when the special
 * form of the decomposition is advantageous.
 *
 * \tparam RandomContainer A random access container
 * \param a \f$N\times N\f$ matrix to be decomposed.  On exit, contains the Householder vectors
 *          \f$u\f$ in its lower triangular part and the strict upper triangular part of \f$R\f$ in
 *          its strict upper triangular part.
 * \param n Rank of the matrix \f$N\f$
 * \param c On exit, contains the normalization \f$\frac{1}{2} u\cdot u\f$ of the Householder
 *          vectors \f$u\f$ .
 * \param d On exit, contains the diagonal part of the \f$R\f$ matrix.
 * \return \c true if the matrix is singular; \c false otherwise.
 *
 * \todo templatize on container element type
 */
template <class RandomContainer>
bool qrdcmp(RandomContainer &a, unsigned n, RandomContainer &c, RandomContainer &d) {
  Require(a.size() == n * n);

  using rtt_dsxx::square;
  using std::sqrt;

  c.resize(n);
  d.resize(n);

  bool singular = false;

  for (unsigned i = 0; i + 1 < n; ++i) {

    // Compute scaling for the ith Householder vector (to prevent overflow)
    double scale = fabs(a[i + n * i]);
    for (unsigned j = i + 1; j < n; j++)
      scale = std::max(scale, fabs(a[j + n * i]));

    if (std::abs(scale) < std::numeric_limits<double>::min()) {

      // ith column is already zeroed from ith element down; the matrix is singular, and the
      // Householder vector is also zero.
      c[i] = d[i] = 0;
      singular = true;
    } else {

      // Compute the Householder vector.
      double sigma = square(a[i + n * i] /= scale);
      for (unsigned j = i + 1; j < n; j++)
        sigma += square(a[j + n * i] /= scale);
      sigma = sqrt(sigma);
      if (a[i + n * i] < 0.0)
        sigma = -sigma;
      // choose sign to minimize roundoff

      // Compute Q*A
      a[i + n * i] += sigma;
      c[i] = sigma * a[i + n * i];
      d[i] = -scale * sigma;
      for (unsigned j = i + 1; j < n; j++) {
        double sum = a[i + n * i] * a[i + n * j];
        for (unsigned k = i + 1; k < n; k++)
          sum += a[k + n * i] * a[k + n * j];
        sum /= -c[i];
        for (unsigned k = i; k < n; k++)
          a[k + n * j] += a[k + n * i] * sum;
      }
    }
  }
  d[n - 1] = a[n - 1 + n * (n - 1)];
  if (std::abs(d[n - 1]) < std::numeric_limits<double>::min() || !rtt_dsxx::isFinite(d[n - 1]))
    singular = true;

  Ensure(a.size() == n * n);
  Ensure(c.size() == n);
  Ensure(d.size() == n);
  return singular;
}

} // end namespace rtt_linear

#endif // linear_qrdcmp_t_hh

//------------------------------------------------------------------------------------------------//
// end of linear/qrdcmp.t.hh
//------------------------------------------------------------------------------------------------//
