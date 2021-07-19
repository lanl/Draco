//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   linear/btridag.t.hh
 * \author Kent Budge
 * \date   Wed Sep 15 13:03:41 MDT 2010
 * \brief  Implementation of block tridiagonal solver
 * \note   Copyright (C) 2016-2021 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef linear_btridag_t_hh
#define linear_btridag_t_hh

#include "ds++/Slice.hh"
#include "linear/btridag.hh"
#include "linear/ludcmp.hh"
#include <vector>

namespace rtt_linear {

//------------------------------------------------------------------------------------------------//
/*!
 * Solve a block tridiagonal system of linear equations.
 *
 * \tparam FieldVector A random access container whose element type must represent a field, such as
 *         the reals or the complex numbers.
 *
 * \param a Subdiagonal of coefficient matrix.
 * \param b Diagonal of coefficient matrix.
 * \param c Superdiagonal of coefficient matrix.
 * \param r Right-hand side of system of equations.
 * \param n Number of blocks
 * \param m Number of variables per block
 * \param u On return, contains the solution of the system.
 *
 * \throw std::range_error If the system is not diagonal dominant.
 */
template <class FieldVector>
void btridag(FieldVector const &a, FieldVector const &b, FieldVector const &c, FieldVector const &r,
             unsigned const n, unsigned const m, FieldVector &u) {
  Require(a.size() == n * m * m);
  Require(b.size() == n * m * m);
  Require(c.size() == n * m * m);
  Require(r.size() == n * m);
  Require(u.size() == n * m);

  using namespace rtt_dsxx;
  using Field = typename FieldVector::value_type;

  if (n == 0) {
    return;
  } else {
    std::vector<Field> gam(n * m * m);
    std::vector<Field> piv(m * m);
    std::vector<Field> rbet(b.begin(), b.begin() + m * m);
    std::vector<unsigned> indx(m);
    double d;
    ludcmp(rbet, indx, d);
    std::copy(r.begin(), r.begin() + m, u.begin());
    Slice<typename FieldVector::iterator> u0(u.begin(), m);
    lubksb(rbet, indx, u0);
    for (unsigned j = 1; j < n; ++j) {
      std::copy(c.begin() + (j - 1) * m * m, c.begin() + j * m * m, gam.begin() + j * m * m);

      for (unsigned k = 0; k < m; ++k) {
        Slice<typename FieldVector::iterator> gamj(gam.begin() + j * m * m + k * m, m);

        lubksb(rbet, indx, gamj);
      }

      for (unsigned k1 = 0; k1 < m; ++k1) {
        Field sumr = r[j * m + k1];
        for (unsigned k2 = 0; k2 < m; ++k2) {
          sumr -= a[j * m * m + k1 + m * k2] * u[(j - 1) * m + k2];
          Field sumb = b[j * m * m + k1 + m * k2];
          for (unsigned k3 = 0; k3 < m; ++k3) {
            sumb -= a[j * m * m + k1 + m * k3] * gam[j * m * m + k3 + m * k2];
          }
          rbet[k1 + m * k2] = sumb;
        }
        u[j * m + k1] = sumr;
      }
      ludcmp(rbet, indx, d);
      Slice<typename FieldVector::iterator> uj(u.begin() + j * m, m);
      lubksb(rbet, indx, uj);
    }
    for (unsigned j = n - 2; j < n - 1; --j) {
      for (unsigned k1 = 0; k1 < m; ++k1) {
        Field sumr = u[j * m + k1];
        for (unsigned k2 = 0; k2 < m; ++k2) {
          sumr -= gam[(j + 1) * m * m + k1 + m * k2] * u[(j + 1) * m + k2];
        }
        u[j * m + k1] = sumr;
      }
    }
  }
}

} // end namespace rtt_linear

#endif // linear_btridag_t_hh

//------------------------------------------------------------------------------------------------//
// end of btridag.t.hh
//------------------------------------------------------------------------------------------------//
