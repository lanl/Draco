//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   linear/tqli.t.hh
 * \author Kent Budge
 * \date   Thu Sep  2 15:00:32 2004
 * \brief  Find eigenvectors and eigenvalues of a symmetric matrix that has been reduced to
 *         tridiagonal form via a call to tred2.
 * \note   Copyright (C) 2004-2021 Triad National Security, LLC.,  All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef linear_tqli_t_hh
#define linear_tqli_t_hh

#include "ds++/Assert.hh"
#include "ds++/DracoMath.hh"
#include "linear/tqli.hh"
#include <stdexcept>

namespace rtt_linear {

//------------------------------------------------------------------------------------------------//
/*!
 * Find eigenvectors and eigenvalues of a symmetric matrix that has been reduced to tridiagonal form
 * via a call to rtt_linear::tred2.
 *
 * \tparam FieldVector1 A random access container on a field type.
 * \tparam FieldVector2 A random access container on a field type.
 * \tparam FieldVector3 A random access container on a field type.
 *
 * \param[in,out] d Diagonal of the matrix.  On return, the eigenvalues of the matrix.
 * \param[in] e Superdiagonal of the matrix.
 * \param[in] n Order of the matrix.
 * \param[in,out] z The rotation matrix (to tridiagonal form) calculated by rtt_linear::tred2.  On
 *                  return, the eigenvectors of the matrix.
 *
 * If the matrix is tridiagonal to begin with, then z should be set to the identity matrix.
 *
 * \pre \c d.size()==n
 * \pre \c e.size()==n
 * \pre \c z.size()==n*n
 */
template <class FieldVector1, class FieldVector2, class FieldVector3>
void tqli(FieldVector1 &d, FieldVector2 &e, const unsigned n, FieldVector3 &z) {
  Require(d.size() == n);
  Require(e.size() == n);
  Require(z.size() == n * n);

  using namespace std;
  using namespace rtt_dsxx;

  // minimum representable value
  double const mrv = std::numeric_limits<typename FieldVector1::value_type>::min();

  for (unsigned i = 1; i < n; ++i) {
    e[i - 1] = e[i];
  }
  e[n - 1] = 0.0;
  for (unsigned l = 0; l < n; ++l) {
    unsigned iter = 0;
    unsigned m;
    do {
      for (m = l; m + 1 < n; ++m) {
        const double dd = abs(d[m]) + abs(d[m + 1]);
        if (rtt_dsxx::soft_equiv(abs(e[m]) + dd, dd))
          break;
      }
      if (m != l) {
        if (iter++ == 300) {
          throw range_error("tqli: no convergence");
        }
        double g = (d[l + 1] - d[l]) / (2.0 * e[l]);
        double r = pythag(g, 1.0);
        g = d[m] - d[l] + e[l] / (g + sign(r, g));
        double s = 1.0;
        double c = 1.0;
        double p = 0.0;
        unsigned i;
        Check(m > 0);
        for (i = m - 1; i >= l && i < m; --i) {
          double f = s * e[i];
          const double b = c * e[i];
          e[i + 1] = (r = pythag(f, g));
          if (std::abs(r) < mrv) {
            d[i + 1] -= p;
            e[m] = 0.0;
            break;
          }
          s = f / r;
          c = g / r;
          g = d[i + 1] - p;
          r = (d[i] - g) * s + 2.0 * c * b;
          p = s * r;
          d[i + 1] = g + p;
          g = c * r - b;
          for (unsigned k = 0; k < n; k++) {
            f = z[k + n * (i + 1)];
            z[k + n * (i + 1)] = s * z[k + n * i] + c * f;
            z[k + n * i] = c * z[k + n * i] - s * f;
          }
        }
        if (std::abs(r) < mrv && i >= l)
          continue;
        d[l] -= p;
        e[l] = g;
        e[m] = 0.0;
      }
    } while (m != l);
  }
}

} // end namespace rtt_linear

#endif // linear_tqli_t_hh

//------------------------------------------------------------------------------------------------//
// end of linear/tqli.t.hh
//------------------------------------------------------------------------------------------------//
