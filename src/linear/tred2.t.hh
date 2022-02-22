//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   linear/tred2.t.hh
 * \author Kent Budge
 * \date   Thu Sep  2 14:49:55 2004
 * \brief  Householder reduce a symmetric matrix.
 * \note   Copyright (C) 2021-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef linear_tred2_t_hh
#define linear_tred2_t_hh

#include "ds++/Assert.hh"
#include "ds++/DracoMath.hh"
#include "linear/tred2.hh"

namespace rtt_linear {

// Use explicit instantiations.

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Householder-reduce a symmetric matrix
 *
 * \tparam FieldVector1 A random access container on a field type.
 * \tparam FieldVector2 A random access container on a field type.
 * \tparam FieldVector3 A random access container on a field type.
 *
 * \param[in,out] a Symmetric matrix stored in [r+n*c] form, that is, as a full matrix.  On return,
 *                  this is replaced by the rotation matrix used to effect the reduction.  This is
 *                  needed for any subsequent call to tqli if the eigenvectors are desired.
 * \param[in] n Dimension of the matrix
 * \param[out] d Diagonal of reduced matrix
 * \param[out] e Superdiagonal of reduced matrix
 *
 * \pre \c a.size()==n*n
 * \pre \c Is_Symmetric(a,n)
 * \pre \c n>0
 *
 * \post \c a.size()==n*n
 * \post \c d.size()==n
 * \post \c e.size()==n
 */
template <typename FieldVector1, typename FieldVector2, typename FieldVector3>
void tred2(FieldVector1 &a, unsigned n, FieldVector2 &d, FieldVector3 &e) {
  Require(a.size() == n * n);
  // O(N*N)    Require(is_symmetric_matrix(a,n));
  Require(n > 0);

  using namespace rtt_dsxx;
  using Field = typename FieldVector1::value_type;

  // minimum representable value
  double const mrv = std::numeric_limits<Field>::min();

  d.resize(n);
  e.resize(n);

  for (unsigned i = n - 1; i > 0; i--) {
    const unsigned l = i - 1;
    Field h = 0.0;
    Field scale = 0.0;
    if (l > 0) {
      for (unsigned k = 0; k <= l; k++) {
        scale += std::abs(a[i + n * k]);
      }
      if (std::abs(scale) < mrv) {
        e[i] = a[i + n * l];
      } else {
        for (unsigned k = 0; k <= l; k++) {
          a[i + n * k] /= scale;
          h += a[i + n * k] * a[i + n * k];
        }
        Field f = a[i + n * l];
        Field g = -sign(std::sqrt(h), f);
        e[i] = scale * g;
        h -= f * g;
        a[i + n * l] = f - g;
        f = 0.0;
        for (unsigned j = 0; j <= l; j++) {
          a[j + n * i] = a[i + n * j] / h;
          g = 0.0;
          for (unsigned k = 0; k <= j; k++) {
            g += a[j + n * k] * a[i + n * k];
          }
          for (unsigned k = j + 1; k <= l; k++) {
            g += a[k + n * j] * a[i + n * k];
          }
          e[j] = g / h;
          f += e[j] * a[i + n * j];
        }
        const double hh = f / (h + h);
        for (unsigned j = 0; j <= l; j++) {
          f = a[i + n * j];
          e[j] = g = e[j] - hh * f;
          for (unsigned k = 0; k <= j; k++) {
            a[j + n * k] -= f * e[k] + g * a[i + n * k];
          }
        }
      }
    } else {
      e[i] = a[i + n * l];
    }
    d[i] = h;
  }
  d[0] = 0.0;
  e[0] = 0.0;
  for (unsigned i = 0; i < n; i++) {
    if (std::abs(d[i]) > mrv) {
      for (unsigned j = 0; j < i; j++) {
        double g = 0.0;
        for (unsigned k = 0; k < i; k++) {
          g += a[i + n * k] * a[k + n * j];
        }
        for (unsigned k = 0; k < i; k++) {
          a[k + n * j] -= g * a[k + n * i];
        }
      }
    }
    d[i] = a[i + n * i];
    a[i + n * i] = 1.0;
    for (unsigned j = 0; j < i; j++) {
      a[j + n * i] = a[i + n * j] = 0.0;
    }
  }
}

} // end namespace rtt_linear

#endif // linear_tred2_t_hh

//------------------------------------------------------------------------------------------------//
// end of linear/tred2.t.hh
//------------------------------------------------------------------------------------------------//
