//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   quadrature/gauleg.hh
 * \author Kent Budge
 * \date   Tue Sep 14 13:16:09 2004
 * \brief  Gauss-Legendre quadrature
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_quadrature_gauleg_hh
#define rtt_quadrature_gauleg_hh

#include "ds++/Soft_Equivalence.hh"
#include "units/PhysicalConstants.hh"
#include <gsl/gsl_sf_legendre.h>
#include <limits>

namespace rtt_quadrature {
//------------------------------------------------------------------------------------------------//
/*!
 * \brief Gauss-Legendre quadrature
 *
 * Calculate abscissae and weights for Gauss-Legendre quadrature:
 * \f[
 * \int_{x_1}^{x_2}{f(x) dx} = \sum\limits_{j=0}^{N - 1}w_j f(x_j)
 * \f]
 *
 * We will use Newton's method of root finding to determine the
 * abscissas. This algorithm is a modified form of that found in "Numerical
 * Recipes in C."
 *
 * The abcissas are the roots of this recurrence relation:
 * \f[
 * (j+1)P_{j+1} = (2j-1)xP_j - jP_{j-1}
 * \f]
 *
 * The weights are determined from the relation:
 * \f[
 * w_j = \frac{2}{(1-x^2_j)[P'_N(x_j)]^2}
 * \f]
 *
 * This routine scales the range of integration from \f$(x_1,x_2)\f$ to (-1,1)
 * and provides the abscissas \f$ x_j\f$ and weights \f$ w_j \f$ for the
 * Gaussian formula provided above.
 *
 * \tparam FieldVector A random access container on a field type.
 *
 * \param x1 Start of integration interval
 * \param x2 End of integration interval
 * \param x On return, contains abscissae \f$x_j\f$ for quadrature.
 * \param w On return, contains weights \f$w_j\f$ for quadrature.
 * \param n Number of points in quadrature.
 */
template <typename FieldVector>
void gauleg(double const x1, // expect FieldVector::value_type to be promoted to double.
            double const x2, // expect FieldVector::value_type to be promoted to double.
            FieldVector &x, FieldVector &w, unsigned const n) {
  using rtt_dsxx::soft_equiv;
  using rtt_units::PI;
  using std::cos;
  using std::numeric_limits;

  using Field = typename FieldVector::value_type;

  Require(n > 0);
  Require(x2 > x1);

  // convergence tolerance
  Field const tolerance(100 * numeric_limits<Field>::epsilon());

  x.resize(n);
  w.resize(n);

  // number of Gauss points in the half range.
  // The roots are symmetric in the interval.  We only need to search for
  // half of them.
  unsigned const numHrGaussPoints((n + 1) / 2);

  // mid-point of integration range
  Field const mu_m(0.5 * (x2 + x1));
  // length of half the integration range.
  Field const mu_l(0.5 * (x2 - x1));

  // Loop over the desired roots.
  for (size_t iroot = 0; iroot < numHrGaussPoints; ++iroot) {
    // Approximate the i-th root.
    Field z(cos(PI * (static_cast<double>(iroot) + 0.75) / (static_cast<double>(n) + 0.5)));

    // Temporary storage;
    Field z1, pp;

    do // Use Newton's method to refine the value for the i-th root.
    {
      // Evaluate the Legendre polynomials evaluated at z.
      Field p1 = gsl_sf_legendre_Pl(n, z);
      Field p2 = gsl_sf_legendre_Pl(n - 1, z);

      // p1 is now the desired Legendre polynomial evaluated at z. We
      // next compute pp, its derivative, by a standard relation
      // involving also p2, the polynomial of one lower order.
      pp = n * (z * p1 - p2) / (z * z - 1.0);

      // Update via Newton's Method
      z1 = z;
      z = z1 - p1 / pp;

    } while (!soft_equiv(z, z1, tolerance));

    // Roots will be between -1 and 1.0 and symmetric about the origin.
    size_t const idxSymPart(n - iroot - 1);

    // Now, scale the root to tthe desired interval and put in its
    // symmetric counterpart.
    x[iroot] = mu_m - mu_l * z;
    x[idxSymPart] = mu_m + mu_l * z;

    // Compute the associated weight and its symmetric counterpart.
    w[iroot] = 2 * mu_l / ((1 - z * z) * pp * pp);
    w[idxSymPart] = w[iroot];
  }
  return;
}

} // end namespace rtt_quadrature

#endif // rtt_quadrature_gauleg_hh

//------------------------------------------------------------------------------------------------//
// end of quadrature/gauleg.hh
//------------------------------------------------------------------------------------------------//
