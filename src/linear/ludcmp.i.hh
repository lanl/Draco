//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   linear/ludcmp.i.hh
 * \author Kent Budge
 * \date   Thu Jul  1 10:54:20 2004
 * \brief  Implementation of methods of ludcmp.hh
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "ludcmp.hh"
#include "ds++/Assert.hh"
#include "ds++/DracoMath.hh"
#include "ds++/Soft_Equivalence.hh"
#include <stdexcept>
#include <vector>

namespace rtt_linear {

//------------------------------------------------------------------------------------------------//
/*!
 * \brief LU-decompose a nonsingular matrix.
 *
 * \arg \a FieldVector1 A random-access container type on a field.
 * \arg \a IntVector A random-access container type on an integral type.
 *
 * \param a Matrix to decompose.  On return, contains the decomposition.
 * \param indx On return, contains the pivoting map.
 * \param d On return, contains the sign of the determinant.
 *
 * \pre \c a.size()==indx.size()*indx.size()
 */
template <typename FieldVector, typename IntVector>
void ludcmp(FieldVector &a, IntVector &indx, typename FieldVector::value_type &d) {
  Require(a.size() == indx.size() * indx.size());

  using Field = typename FieldVector::value_type;

  Check(indx.size() < UINT_MAX);
  auto const n = static_cast<unsigned>(indx.size());

  std::vector<Field> vv(n);

  d = 1.0;
  for (unsigned i = 0; i < n; ++i) {
    Field big = 0.0;
    for (unsigned j = 0; j < n; ++j) {
      Field const temp = std::abs(a[i + n * j]);
      if (temp > big) {
        big = temp;
      }
    }
    vv[i] = 1.0 / big;
    if (std::isinf(vv[i])) {
      throw std::domain_error("ludcmp:  singular matrix");
    }
  }
  for (unsigned j = 0; j < n; ++j) {
    for (unsigned i = 0; i < j; ++i) {
      Field sum = a[i + n * j];
      for (unsigned k = 0; k < i; ++k) {
        sum -= a[i + n * k] * a[k + n * j];
      }
      a[i + n * j] = sum;
    }
    Field big = 0.0;
    unsigned imax(0);
    for (unsigned i = j; i < n; ++i) {
      Field sum = a[i + n * j];
      for (unsigned k = 0; k < j; ++k) {
        sum -= a[i + n * k] * a[k + n * j];
      }
      a[i + n * j] = sum;
      Field const dum = vv[i] * std::abs(sum);
      if (dum >= big) {
        big = dum;
        imax = i;
      }
    }
    if (j != imax) {
      for (unsigned k = 0; k < n; ++k) {
        Field dum = a[imax + n * k];
        a[imax + n * k] = a[j + n * k];
        a[j + n * k] = dum;
      }
      d = -d;
      vv[imax] = vv[j];
    }
    indx[j] = imax;
    if (j != n - 1) {
      Field dum = 1.0 / a[j + n * j];
      if (std::isinf(dum)) {
        throw std::domain_error("ludcmp:  singular matrix");
      }
      for (unsigned i = j + 1; i < n; ++i)
        a[i + n * j] *= dum;
    }
  }
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Solve the system \f$Ax=b\f$
 *
 * \arg \a FieldVector1 A random-access container type on a field.
 * \arg \a IntVector A random-access container type on an integral type.
 * \arg \a FieldVector2 A random-access container type on a field.
 *
 * \param a LU decomposition of \f$A\f$.
 * \param indx Pivot map for decomposition of \f$A\f$.
 * \param b Right-hand side \f$b\f$.  On return, contains solution \f$x\f$.
 *
 * \pre \c a.size()==indx.size()*indx.size()
 * \pre \c b.size()==indx.size()
 */
template <typename FieldVector1, typename IntVector, typename FieldVector2>
void lubksb(FieldVector1 const &a, IntVector const &indx, FieldVector2 &b) {
  Require(a.size() == indx.size() * indx.size());
  Require(b.size() == indx.size());

  using Field = typename FieldVector2::value_type;

  // minimum representable value
  double const mrv = std::numeric_limits<Field>::min();
  Check(indx.size() < UINT_MAX);
  auto const n = static_cast<unsigned>(indx.size());

  unsigned ii = 0;

  for (unsigned i = 0; i < n; ++i) {
    unsigned ip = indx[i];
    Field sum = b[ip];
    b[ip] = b[i];
    if (ii != 0) {
      for (unsigned j = ii - 1; j < i; ++j)
        sum -= a[i + n * j] * b[j];
    } else {
      if (std::abs(sum) > mrv)
        ii = i + 1;
    }
    b[i] = sum;
  }
  for (unsigned i = n - 1; i < n; --i) {
    Field sum = b[i];
    for (unsigned j = i + 1; j < n; ++j)
      sum -= a[i + n * j] * b[j];
    b[i] = sum / a[i + n * i];
  }
}

} // end namespace rtt_linear

//------------------------------------------------------------------------------------------------//
// end of ludcmp.i.hh
//------------------------------------------------------------------------------------------------//
