//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   linear/qrupdt.t.hh
 * \author Kent Budge
 * \date   Tue Aug 10 11:59:48 2004
 * \brief  Update the QR decomposition of a square matrix
 * \note   Copyright (C) 2016-2021 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef linear_qrupdt_t_hh
#define linear_qrupdt_t_hh

#include "qrupdt.hh"
#include "rotate.hh"
#include "ds++/Assert.hh"
#include "ds++/DracoMath.hh"

namespace rtt_linear {
//------------------------------------------------------------------------------------------------//
/*!
 * \brief Update the QR decomposition of a square matrix.
 *
 * One of the great advantages of the QR decomposition is the ease with which it can be updated.  If
 * \f$A=QR\f$ where \f$Q\f$ is an orthonormal matrix and \f$R\f$ is upper triangular, then
 * \f$A+s\bigotimes t = Q(R+u \bigotimes v)\f$ where \f$v=t\f$ and \f$u=Q^Ts\f$.  Updating the QR
 * decomposition takes of order \f$N^2\f$ operations rather than the \f$N^3\f$ operations of a full
 * matrix inversion.
 *
 * \tparam RandomContainer A random access container.
 *
 * \param r Upper triangular matrix of the QR decomposition.
 * \param qt Transpose of orthonormal matrix of the QR decomposition.
 * \param n Rank of the matrix.
 * \param u Update vector \f$u\f$.  Destroyed on return.
 * \param v Update vector \f$v\f$.  Destroyed on return.
 *
 * \todo Templatize on container element type
 */
template <class RandomContainer>
void qrupdt(RandomContainer &r, RandomContainer &qt, const unsigned n, RandomContainer &u,
            RandomContainer &v) {

  Require(r.size() == n * n);
  Require(qt.size() == n * n);
  Require(u.size() == n);
  Require(v.size() == n);

  using std::fabs;
  using namespace rtt_dsxx;

  // minumum representable value
  double const mrv = std::numeric_limits<typename RandomContainer::value_type>::min();

  // Find first nonzero element of u.
  int k;
  for (k = n - 1; k >= 0; --k) {
    if (std::abs(u[k]) > mrv)
      break;
  }
  if (k < 0)
    k = 0;
  for (int i = k - 1; i >= 0; i--) {
    rotate(r, qt, n, i, u[i], -u[i + 1]);
    if (std::abs(u[i]) < mrv) {
      u[i] = fabs(u[i + 1]);
    } else if (fabs(u[i]) > fabs(u[i + 1])) {
      u[i] = fabs(u[i]) * sqrt(1.0 + square(u[i + 1] / u[i]));
    } else {
      u[i] = fabs(u[i + 1]) * sqrt(1.0 + square(u[i] / u[i + 1]));
    }
  }
  for (unsigned j = 0; j < n; j++)
    r[0 + n * j] += u[0] * v[j];
  for (int i = 0; i < k; i++)
    rotate(r, qt, n, i, r[i + n * i], -r[i + 1 + n * i]);

  Ensure(r.size() == n * n);
  Ensure(qt.size() == n * n);
}

} // end namespace rtt_linear

#endif // linear_qrupdt_t_hh

//------------------------------------------------------------------------------------------------//
// end of linear/qrupdt.t.hh
//------------------------------------------------------------------------------------------------//
