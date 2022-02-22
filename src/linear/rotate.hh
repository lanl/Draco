//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   linear/rotate.hh
 * \author Kent Budge
 * \date   Tue Aug 10 12:37:43 2004
 * \brief  Perform a Jacobi rotation on a pair of matrices.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef linear_rotate_hh
#define linear_rotate_hh

#include "ds++/Assert.hh"
#include "ds++/Soft_Equivalence.hh"
#include <cmath>
#include <limits>

namespace rtt_linear {

using std::sqrt;

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Perform a Jacobi rotation on a QR decomposition.
 *
 * Perform a Jacobi rotation on rows \f$i\f$ and \f$i+1\f$ of the orthonormal matrix \f$Q\f$ and the
 * upper triangular matrix \f$R\f$.  The rotation is parametrized by \f$a\f$ and \f$b\f$ which
 * satisfy \f$\cos\theta = \frac{a}{\sqrt{a^2+b^2}}\f$ and \f$\sin\theta =
 * \frac{a}{\sqrt{a^2+b^2}}\f$
 *
 * \arg \a RandomContainer A random access container.
 *
 * \param r Upper triangular matrix on which to perform Jacobi rotation.
 * \param qt Orthonormal matrix on which to perform Jacobi rotations.
 * \param n Rank of the matrix
 * \param i Row of the matrix on which to perform the Jacobi rotation.
 * \param a First rotation parameter
 * \param b Second rotation parameter
 *
 * \pre \c r.size()==n*n
 * \pre \c qt.size()==n*n
 * \pre \c i+1<n
 *
 * \post \c r.size()==n*n
 * \post \c qt.size()==n*n
 */
template <class RandomContainer>
void rotate(RandomContainer &r, RandomContainer &qt, const unsigned n, unsigned i, double a,
            double b) {
  Require(r.size() == n * n);
  Require(qt.size() == n * n);
  Require(i + 1 < n);

  using std::fabs;

  // cosine and sine of rotation
  double c, s;
  // minimum representable value
  double const mrv = std::numeric_limits<typename RandomContainer::value_type>::min();

  if (std::abs(a) < mrv) {
    c = 0.0;
    s = (b > 0.0 ? 1.0 : -1.0);
  } else if (fabs(a) > fabs(b)) {
    double fact = b / a;
    c = 1.0 / sqrt(1.0 + fact * fact);
    if (a < 0)
      c = -c;
    s = fact * c;
  } else {
    double fact = a / b;
    s = 1.0 / sqrt(1.0 + fact * fact);
    if (b < 0)
      s = -s;
    c = fact * s;
  }
  for (unsigned j = i; j < n; j++) {
    double y = r[i + n * j];
    double w = r[i + 1 + n * j];
    r[i + n * j] = c * y - s * w;
    r[i + 1 + n * j] = s * y + c * w;
  }
  for (unsigned j = 0; j < n; j++) {
    double y = qt[i + n * j];
    double w = qt[i + 1 + n * j];
    qt[i + n * j] = c * y - s * w;
    qt[i + 1 + n * j] = s * y + c * w;
  }

  Ensure(r.size() == n * n);
  Ensure(qt.size() == n * n);
}

} // end namespace rtt_linear

#endif // linear_rotate_hh

//------------------------------------------------------------------------------------------------//
// end of linear/rotate.hh
//------------------------------------------------------------------------------------------------//
