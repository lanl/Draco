//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   linear/svdcmp.t.hh
 * \author Kent Budge
 * \date   Mon Aug  9 13:17:31 2004
 * \brief  Calculate the singular value decomposition of a matrix.
 * \note   Copyright (C) 2021-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef linear_svdcmp_t_hh
#define linear_svdcmp_t_hh

#include "ds++/DracoMath.hh"
#include "linear/svdcmp.hh"
#include <sstream>
#include <vector>

namespace rtt_linear {
//------------------------------------------------------------------------------------------------//
/*!
 * \brief Compute the singular value decomposition of a matrix.
 *
 * Compute the decomposition of a matrix \f$ A=UWV^T \f$ where \f$ U \f$ has the same shape as the
 * original matrix; \f$ W \f$ is diagonal with rank equal to the column order of \f$ A \f$; and \f$
 * V \f$ is a square full matrix of rank equal to the column order of \f$ A \f$.
 *
 * The singular value decomposition is tremendously useful for manipulation both nonsquare matrices
 * and nearly singular square matrices.  The following routine is very robust.
 *
 * \tparam RandomContainer A random access container type
 * 
 * \param a Matrix to be decomposed.  On exit, contains \f$ U \f$.
 * \param m Number of rows in a
 * \param n Number of columns in a
 * \param w On exit, contains \f$ W \f$.
 * \param v On exit, contains \f$ V \f$.
 * \todo Templatize on container element type
 */
template <typename RandomContainer>
void svdcmp(RandomContainer &a, const unsigned m, const unsigned n, RandomContainer &w,
            RandomContainer &v) {
  Require(a.size() == m * n);

  using namespace rtt_dsxx;
  using std::fabs;
  using std::max;
  using std::min;
  using std::sqrt;

  // More than 30 iterations says something is terribly wrong -- this shouldn't happen even for very
  // large matrices.
  const unsigned MAX_ITERATIONS = 30;
  // minimum representable value
  double const mrv = std::numeric_limits<typename RandomContainer::value_type>::min();
  double const eps = std::numeric_limits<typename RandomContainer::value_type>::epsilon();

  w.resize(n);
  v.resize(n * n);
  std::vector<double> rv1(n);

  // Reduce to bidiagonal form.
  double g = 0;
  double scale = 0;
  double norm = 0;
  unsigned l = 0; // in case n==0
  for (unsigned i = 0; i < n; i++) {
    l = i + 1;
    rv1[i] = scale * g;
    g = 0;
    double s = 0.0;
    scale = 0;
    if (i < m) {
      for (unsigned k = i; k < m; k++)
        scale += fabs(a[k + m * i]);
      if (std::abs(scale) > mrv) {
        double rscale = 1 / scale;
        for (unsigned k = i; k < m; k++) {
          a[k + m * i] *= rscale;
          s += square(a[k + m * i]);
        }
        double f = a[i + m * i];
        g = sqrt(s);
        if (f > 0.0)
          g = -g;
        double h = f * g - s;
        a[i + m * i] = f - g;
        for (unsigned j = l; j < n; j++) {
          s = 0;
          for (unsigned k = i; k < m; k++)
            s += a[k + m * i] * a[k + m * j];
          Check(std::abs(h) > mrv);
          f = s / h;
          for (unsigned k = i; k < m; k++)
            a[k + m * j] += f * a[k + m * i];
        }
        for (unsigned k = i; k < m; k++)
          a[k + m * i] *= scale;
      }
    }
    w[i] = scale * g;
    g = 0;
    s = 0;
    scale = 0;
    if (i < m) {
      for (unsigned k = l; k < n; k++)
        scale += fabs(a[i + m * k]);
      if (std::abs(scale) > mrv) {
        double rscale = 1 / scale;
        for (unsigned k = l; k < n; k++) {
          a[i + m * k] *= rscale;
          s += square(a[i + m * k]);
        }
        double f = a[i + m * l];
        g = sqrt(s);
        if (f > 0.0)
          g = -g;
        double h = f * g - s;
        a[i + m * l] = f - g;
        double rh = 1 / h;
        for (unsigned k = l; k < n; k++)
          rv1[k] = a[i + m * k] * rh;
        for (unsigned j = l; j < m; j++) {
          s = 0;
          for (unsigned k = l; k < n; k++)
            s += a[j + m * k] * a[i + m * k];
          for (unsigned k = l; k < n; k++)
            a[j + m * k] += s * rv1[k];
        }
        for (unsigned k = l; k < n; k++)
          a[i + m * k] *= scale;
      }
    }
    norm = max(norm, fabs(w[i]) + fabs(rv1[i]));
  }
  // Accumulation of right-hand transformations
  for (unsigned i = n - 1; i < n; i--) {
    if (i != n - 1) {
      if (std::abs(g) > mrv) {
        double rg = 1 / g;
        for (unsigned j = l; j < n; j++)
          v[j + n * i] = rg * (a[i + m * j] / a[i + m * l]);
        // ordering of above expression important to prevent underflow
        for (unsigned j = l; j < n; j++) {
          double s = 0;
          for (unsigned k = l; k < n; k++)
            s += a[i + m * k] * v[k + n * j];
          for (unsigned k = l; k < n; k++)
            v[k + n * j] += s * v[k + n * i];
        }
      }
      for (unsigned j = l; j < n; j++)
        v[i + n * j] = v[j + n * i] = 0;
    }
    v[i + n * i] = 1;
    g = rv1[i];
    l = i;
  }
  // Accumulation of left-hand transformations
  for (unsigned i = min(m, n) - 1; i < min(m, n); i--) {
    l = i + 1;
    g = w[i];
    for (unsigned j = l; j < n; j++)
      a[i + m * j] = 0.0;
    if (std::abs(g) > mrv) {
      double rg = 1 / g;
      for (unsigned j = l; j < n; j++) {
        double s = 0;
        for (unsigned k = l; k < m; k++)
          s += a[k + m * i] * a[k + m * j];
        double f = rg * (s / a[i + m * i]);
        // ordering of above expression important to prevent underflow
        for (unsigned k = i; k < m; k++)
          a[k + m * j] += f * a[k + m * i];
      }
      for (unsigned j = i; j < m; j++)
        a[j + m * i] *= rg;
    } else {
      for (unsigned j = i; j < m; j++)
        a[j + m * i] = 0;
    }
    a[i + m * i] += 1;
  }

  // Reduce to diagonal form.
  for (unsigned k = n - 1; k < n; k--) {
    unsigned k1 = k - 1;
    unsigned its = 0;
    for (; its < 30; its++) // Allow up to 30 iterations.
    {
      bool flag = true;
      // Check for splitting.
      Check(rtt_dsxx::soft_equiv(rv1[0] + norm, norm));
      unsigned ell = k;
      for (; ell <= k; ell--) {
        unsigned l1 = ell - 1;
        if (rtt_dsxx::soft_equiv(fabs(rv1[ell]) + norm, norm, eps)) {
          flag = false;
          break;
        }
        if (rtt_dsxx::soft_equiv(fabs(w[l1]) + norm, norm, eps))
          break;
      }
      if (flag) {
        unsigned l1 = ell - 1;
        double c = 0;
        double s = 1;
        for (unsigned i = ell; i <= k; i++) {
          double f = s * rv1[i];
          rv1[i] *= c;
          if (!rtt_dsxx::soft_equiv(fabs(f) + norm, norm, eps)) {
            g = w[i];
            double h = pythag(f, g);
            w[i] = h;
            c = g / h;
            s = -f / h;
            for (unsigned j = 0; j < m; j++) {
              double y = a[j + m * l1];
              double z = a[j + m * i];
              a[j + m * l1] = y * c + z * s;
              a[j + m * i] = -y * s + z * c;
            }
          }
        }
      }
      double z = w[k];
      if (ell == k) {
        if (z < 0.0) {
          w[k] = -z;
          for (unsigned j = 0; j < n; j++)
            v[j + n * k] = -v[j + n * k];
        }
        break;
      } else {
        double x = w[ell];
        double y = w[k1];
        g = rv1[k1];
        double h = rv1[k];
        double f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2 * h * y);
        g = pythag(f, 1.0);
        if (f < 0.0)
          g = -g;
        f = ((x - z) * (x + z) + h * (y / (f + g) - h)) / x;
        double c = 1;
        double s = 1;
        for (unsigned j = ell; j < k; j++) {
          int i = j + 1;
          g = rv1[i];
          y = w[i];
          h = s * g;
          g = c * g;
          z = pythag(f, h);
          rv1[j] = z;
          c = f / z;
          s = h / z;
          f = x * c + g * s;
          g = -x * s + g * c;
          h = y * s;
          y *= c;
          for (unsigned jj = 0; jj < n; jj++) {
            x = v[jj + n * j];
            z = v[jj + n * i];
            v[jj + n * j] = x * c + z * s;
            v[jj + n * i] = -x * s + z * c;
          }
          z = pythag(f, h);
          w[j] = z;
          if (std::abs(z) > eps) {
            c = f / z;
            s = h / z;
          }
          f = c * g + s * y;
          x = -s * g + c * y;
          for (unsigned jj = 0; jj < m; jj++) {
            y = a[jj + m * j];
            z = a[jj + m * i];
            a[jj + m * j] = y * c + z * s;
            a[jj + m * i] = -y * s + z * c;
          }
        }
        rv1[ell] = 0;
        rv1[k] = f;
        w[k] = x;
      }
    }
    if (its == MAX_ITERATIONS) {
      std::ostringstream message;
      message << "svdcmp: no convergence for singular value " << k;
      throw std::range_error(message.str());
    }
  }
}

} // end namespace rtt_linear

#endif // linear_svdcmp_t_hh

//------------------------------------------------------------------------------------------------//
// end of linear/svdcmp.t.hh
//------------------------------------------------------------------------------------------------//
