//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   linear/gaussj.t.hh
 * \author Kent Budge
 * \brief  Solve a linear system by Gaussian elimination.
 * \note   Copyright (C) 2021-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef linear_gaussj_t_hh
#define linear_gaussj_t_hh

#include "ds++/Assert.hh"
#include "ds++/DracoMath.hh"
#include "ds++/Field_Traits.hh"
#include "linear/gaussj.hh"
#include <algorithm>
#include <sstream>
#include <vector>

namespace rtt_linear {
using std::fabs;
using std::max;
using std::sqrt;

//------------------------------------------------------------------------------------------------//
//! Is a double-subscript random container square?
template <class DoubleRandomContainer> bool is_square(DoubleRandomContainer const &A) {
  Check(A.size() < UINT_MAX);
  auto const n = static_cast<unsigned>(A.size());
  for (unsigned i = 0; i < n; ++i) {
    if (A[i].size() != n)
      return false;
  }
  return true;
}

//------------------------------------------------------------------------------------------------//
/*!
 * \tparam RandomContainer A random access container type
 *
 * \param A Coefficient matrix of the system of equations. Destroyed on return.
 * \param n Rank of matrix A
 * \param b Right hand side of the system of equations. Replaced by the solution on return.
 * \param m Column count of the right hand side of the system of equations.  Setting this to a value
 *          other than one amounts to simultaneously solving m systems of equations.
 */
template <class RandomContainer>
void gaussj(RandomContainer &A, unsigned const n, RandomContainer &b, unsigned const m) {
  using namespace std;
  using namespace rtt_dsxx;

  Require(A.size() == n * n);
  Require(b.size() == n * m);
  double const eps = std::numeric_limits<typename RandomContainer::value_type>::epsilon();

  vector<int> indxc(n);
  vector<int> indxr(n);
  vector<int> ipiv(n, 0);

  unsigned irow(0), icol(0);
  for (unsigned i = 0; i < n; i++) {
    double big = 0.0;
    for (unsigned j = 0; j < n; j++) {
      if (ipiv[j] != 1) {
        for (unsigned k = 0; k < n; k++) {
          if (ipiv[k] == 0) {
            if (fabs(value(A[j + n * k])) >= big) {
              big = fabs(value(A[j + n * k]));
              irow = j;
              icol = k;
            }
          }
        }
      }
    }
    ++ipiv[icol];
    if (irow != icol) {
      for (unsigned l = 0; l < n; l++) {
        swap(A[irow + n * l], A[icol + n * l]);
      }
      for (unsigned l = 0; l < m; l++) {
        swap(b[irow + n * l], b[icol + n * l]);
      }
    }
    indxr[i] = irow;
    indxc[i] = icol;
    if (rtt_dsxx::soft_equiv(A[icol + n * icol], 0.0, eps)) {
      throw invalid_argument("gaussj:  singular matrix");
    }
    double const pivinv = 1.0 / A[icol + n * icol];
    A[icol + n * icol] = 1.0;
    for (unsigned l = 0; l < n; l++) {
      A[icol + n * l] *= pivinv;
    }
    for (unsigned l = 0; l < m; l++) {
      b[icol + n * l] *= pivinv;
    }
    for (unsigned ll = 0; ll < n; ll++) {
      if (ll != icol) {
        double const dum = A[ll + n * icol];
        A[ll + n * icol] = 0.0;
        for (unsigned l = 0; l < n; l++) {
          A[ll + n * l] -= A[icol + n * l] * dum;
        }
        for (unsigned l = 0; l < m; l++) {
          b[ll + n * l] -= b[icol + n * l] * dum;
        }
      }
    }
  }
  for (unsigned l = n - 1; l < n; l--) {
    if (indxr[l] != indxc[l]) {
      for (unsigned k = 0; k < n; k++) {
        swap(A[k + n * indxr[l]], A[k + n * indxc[l]]);
      }
    }
  }

  Ensure(b.size() == n * m);
}

//------------------------------------------------------------------------------------------------//
/*!
 * \tparam DoubleRandomContainer A double-subscript random access container type
 * \tparam RandomContainer A random access container type
 *
 * \param A Coefficient matrix of the system of equations. Destroyed on return.
 * \param b Right hand side of the system of equations. Replacec by the solution of the system on
 *          return.
 */
template <class DoubleRandomContainer, class RandomContainer>
void gaussj(DoubleRandomContainer &A, RandomContainer &b) {
  using namespace std;
  using namespace rtt_dsxx;

  Require(is_square(A));
  Require(b.size() == 0 || b.size() == A.size());

  // minimum representable value
  double const mrv = std::numeric_limits<typename RandomContainer::value_type>::min();
  Check(A.size() < UINT_MAX);
  auto const n = static_cast<unsigned>(A.size());

  vector<int> indxc(n);
  vector<int> indxr(n);
  vector<int> ipiv(n, 0);

  unsigned irow(0), icol(0);
  for (unsigned i = 0; i < n; i++) {
    double big = 0.0;
    for (unsigned j = 0; j < n; j++) {
      if (ipiv[j] != 1) {
        for (unsigned k = 0; k < n; k++) {
          if (ipiv[k] == 0) {
            if (fabs(value(A[j][k])) >= big) {
              big = fabs(value(A[j][k]));
              irow = j;
              icol = k;
            }
          }
        }
      }
    }
    if (fabs(big) < mrv) {
      throw invalid_argument("gaussj:  singular matrix");
    }
    ++ipiv[icol];
    if (irow != icol) {
      for (unsigned l = 0; l < n; l++) {
        swap(A[irow][l], A[icol][l]);
      }
      swap(b[irow], b[icol]);
    }
    indxr[i] = irow;
    indxc[i] = icol;
    if (fabs(value(A[icol][icol])) < mrv) {
      throw invalid_argument("gaussj:  singular matrix");
    }
    double const pivinv = 1.0 / A[icol][icol];
    A[icol][icol] = 1.0;
    for (unsigned l = 0; l < n; l++) {
      A[icol][l] *= pivinv;
    }
    b[icol] *= pivinv;
    for (unsigned ll = 0; ll < n; ll++) {
      if (ll != icol) {
        double const dum = A[ll][icol];
        A[ll][icol] = 0.0;
        for (unsigned l = 0; l < n; l++) {
          A[ll][l] -= A[icol][l] * dum;
        }
        b[ll] -= b[icol] * dum;
      }
    }
  }
  for (unsigned l = n - 1; l < n; l--) {
    if (indxr[l] != indxc[l]) {
      for (unsigned k = 0; k < n; k++) {
        swap(A[k][indxr[l]], A[k][indxc[l]]);
      }
    }
  }

  Ensure(b.size() == A.size());
}

} // end namespace rtt_linear

#endif // linear_gaussj_t_hh

//------------------------------------------------------------------------------------------------//
// end of linear/gaussj.t.hh
//------------------------------------------------------------------------------------------------//
