//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   linear/qr_unpack.t.hh
 * \author Kent Budge
 * \date   Mon Aug  9 13:17:31 2004
 * \brief  Compute an explicit representation of a packed QR decomposition.
 * \note   Copyright (C) 2021-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef linear_qr_unpack_t_hh
#define linear_qr_unpack_t_hh

#include "qr_unpack.hh"
#include "ds++/Assert.hh"
#include "ds++/DracoMath.hh"
#include <algorithm>
#include <limits>
#include <sstream>

namespace rtt_linear {
//------------------------------------------------------------------------------------------------//
/*!
 * \brief Compute an explicit representation of a packed QR decomposition.
 *
 * This procedure takes a packed representation of a QR decomposition, as returned by qrdcmp, and
 * computes explicit \f$Q^T\f$ and \f$R\f$ matrices.  \f$Q^T\f$ rather than \f$Q\f$ is returned
 * because this is what is usually required for further computation (as with qrupdt).
 *
 * \tparam RandomContainer A random access container
 * 
 * \param r \f$N\times N\f$ matrix containing, in its lower triangular part, the Householder vectors
 *          \f$u\f$ of the QR decomposition, and, in its strict upper triangular part, the strict
 *          upper triangular part of \f$R\f$.  On exit, contains \f$R\f$.
 * \param n Rank of matrix \f$N\f$
 * \param c Householder vector normalizations \f$\frac{1}{2}u\cdot u\f$.
 * \param d Diagonal elements of \f$R\f$.
 * \param qt On exit, contains \f$Q^T\f$.
 *
 * \pre \c n*n==r.size()
 * \pre \c n==c.size()
 * \pre \c n=d.size()
 * \pre The elements of <code>c</code> must be set to the correct normalization values for the
 *      \f$ u \f$ contained in \c r.
 *
 * \post \c r[i][j]==0 for all \c j<i
 * \post \c r.size()==n*n
 * \post \c qt.size()==r.size()
 */
template <class RandomContainer>
void qr_unpack(RandomContainer &r, const unsigned n, const RandomContainer &c,
               const RandomContainer &d, RandomContainer &qt) {
  Require(r.size() == n * n);
  Require(n == c.size());
  Require(n == d.size());

  // minimum representable value
  double const mrv = std::numeric_limits<typename RandomContainer::value_type>::min();
  qt.resize(n * n);

  for (unsigned i = 0; i < n; ++i) {
    for (unsigned j = 0; j < n; ++j) {
      qt[i + n * j] = 0;
    }
    qt[i + n * i] = 1;
  }
  // Explicitly form Q transpose
  for (unsigned i = 0; i + 1 < n; ++i) {
    if (fabs(c[i]) > mrv) {
      double rscale = -1 / c[i];
      for (unsigned j = 0; j < n; ++j) {
        double sum = r[i + n * i] * qt[i + n * j];
        for (unsigned k = i + 1; k < n; ++k) {
          sum += r[k + n * i] * qt[k + n * j];
        }
        sum *= rscale;
        for (unsigned k = i; k < n; k++) {
          qt[k + n * j] += r[k + n * i] * sum;
        }
      }
    }
  }
  // Explicitly form r
  for (unsigned i = 0; i < n; i++) {
    r[i + n * i] = d[i];
    for (unsigned k = 0; k < i; k++) {
      r[i + n * k] = 0;
    }
  }

  Ensure(r.size() == n * n);
  Ensure(qt.size() == r.size());
}

} // end namespace rtt_linear

#endif // linear_qr_unpack_t_hh

//------------------------------------------------------------------------------------------------//
// end of linear/qr_unpack.t.hh
//------------------------------------------------------------------------------------------------//
