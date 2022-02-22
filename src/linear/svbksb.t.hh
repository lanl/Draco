//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   linear/svbksb.t.hh
 * \author Kent Budge
 * \date   Tue Aug 10 13:08:03 2004
 * \brief  Solve a linear system from its singular value decomposition.
 * \note   Copyright (C) 2021-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef linear_svbksb_t_hh
#define linear_svbksb_t_hh

#include "ds++/Assert.hh"
#include "ds++/Soft_Equivalence.hh"
#include "linear/svbksb.hh"
#include <limits>
#include <vector>

namespace rtt_linear {
//------------------------------------------------------------------------------------------------//
/*!
 * Solve a linear system given its singular value decomposition.
 *
 * Solves the system \f$ UWV^Tx=b \f$ where \f$ U \f$ is an \f$ M\times N \f$ column-orthogonal
 * matrix, \f$ W \f$ is an \f$ N\times N \f$ diagonal matrix, and \f$ V \f$ is an \f$ N\times N \f$
 * orthonormal matrix.  These matrices define the singular value decomposition of a general \f$
 * M\times N \f$ matrix.  \f$ b \f$ is a vector of length \f$ M \f$ containing the right-hand side
 * and \f$ x \f$ is a vector of length \f$ N \f$ into which the solution will be stored.
 *
 * \tparam RandomContainer A random access container
 * \param u Injection matrix of singular-value decomposition of the linear system.
 * \param w Singular value matrix of singular-value decomposition of the linear system.
 * \param v Projection matrix of singular-value decomposition of the linear system.
 * \param m Rows of the original coefficient matrix
 * \param n Columns of the original coefficient matrix
 * \param b Right-hand side of the linear system.
 * \param x On return, contains the solution of the system.
 *
 * \pre \c u.size()==m*n
 * \pre \c w.size()==n
 * \pre \c b.size()==m
 * \pre \c v.size()==n*n
 *
 * \post \c x.size()==n
 * \post \c x satisfies \f$UWVx=b\f$
 */
template <typename RandomContainer>
void svbksb(const RandomContainer &u, const RandomContainer &w, const RandomContainer &v,
            const unsigned m, const unsigned n, const RandomContainer &b, RandomContainer &x) {
  Require(u.size() == m * n);
  Require(w.size() == n);
  Require(b.size() == m);
  Require(v.size() == n * n);

  using value_type = typename RandomContainer::value_type;
  // minimum representable value
  double const mrv = std::numeric_limits<value_type>::min();
  std::vector<value_type> tmp(n);

  for (unsigned i = 0; i < n; i++) {
    if (std::abs(w[i]) > mrv)
    // Exclude singular values.  This is most of the "magic" of singular value decomposition.
    {
      // Multiply the RHS by transpose of U == inverse of U.
      value_type sum = 0.0;
      for (unsigned j = 0; j < m; j++) {
        sum += u[j + m * i] * b[j];
      }
      // Divide by w.
      tmp[i] = sum / w[i];
    } else {
      tmp[i] = 0.0;
    }
  }
  // Now multiply by V == inverse of transpose of V
  x.resize(n);
  for (unsigned i = 0; i < n; i++) {
    value_type sum = 0.0;
    for (unsigned j = 0; j < n; j++) {
      sum += v[i + n * j] * tmp[j];
    }
    x[i] = sum;
  }

  Ensure(x.size() == n);
}

} // end namespace rtt_linear

#endif // linear_svbksb_t_hh

//------------------------------------------------------------------------------------------------//
// end of linear/svbksb.t.hh
//------------------------------------------------------------------------------------------------//
