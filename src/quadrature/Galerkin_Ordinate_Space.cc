//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   quadrature/Galerkin_Ordinate_Space.cc
 * \author Kent Budge
 * \date   Mon Mar 26 16:11:19 2007
 * \brief  Define methods of class Galerkin_Ordinate_Space
 * \note   Copyright (C) 2012-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "Galerkin_Ordinate_Space.hh"
#include "special_functions/Ylm.hh"
#include "units/PhysicalConstants.hh"
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <iomanip>
#include <iostream>

using namespace rtt_units;
using std::vector;

namespace rtt_quadrature {

//------------------------------------------------------------------------------------------------//
vector<Moment> Galerkin_Ordinate_Space::compute_n2lk_1D_(Quadrature_Class /*unused*/,
                                                         unsigned const N) {
  vector<Moment> result;

  // Choose: l= 0, ..., N-1, k = 0
  int constexpr k(0); // k is always zero for 1D.

  for (unsigned ell = 0; ell < N; ++ell)
    result.emplace_back(Moment(ell, k));

  return result;
}

//------------------------------------------------------------------------------------------------//
vector<Moment> Galerkin_Ordinate_Space::compute_n2lk_1Da_(Quadrature_Class /*unused*/,
                                                          unsigned const N) {
  std::vector<Moment> result;

  // Choose: l= 0, ..., N, k = 0, ..., l to eliminate moments even in xi

  for (int ell = 0; ell < static_cast<int>(N); ++ell)
    for (int k = 0; k <= ell; ++k)
      if ((ell + k) % 2 == 0)
        result.emplace_back(Moment(ell, k)); // Eliminate moments even in eta

  return result;
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Creates a mapping between moment index n and the index pair (k,l).
 */
vector<Moment> Galerkin_Ordinate_Space::compute_n2lk_2D_(Quadrature_Class /*unused*/,
                                                         unsigned const N) {
  std::vector<Moment> result;

  // X-Y symmetry

  for (int ell = 0; ell <= static_cast<int>(N); ++ell)
    for (int k = -ell; k <= ell; ++k)
      if (((ell + abs(k)) % 2 == 0) && ((ell < static_cast<int>(N)) || (k < 0)))
        result.emplace_back(Moment(ell, k));

  return result;
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Creates a mapping between moment index n and the index pair (k,l).
 */
vector<Moment> Galerkin_Ordinate_Space::compute_n2lk_2Da_(Quadrature_Class /*unused*/,
                                                          unsigned const N) {
  std::vector<Moment> result;

  // R-Z symmetry
  // Choose: l= 0, ..., N-1, k = 0, ..., l

  for (int ell = 0; ell < static_cast<int>(N); ++ell)
    for (int k = 0; k <= ell; ++k)
      result.emplace_back(Moment(ell, k));

  // Add N and k>0, k odd
  int ell(N);
  for (int k = 1; k <= ell; k += 2)
    result.emplace_back(Moment(ell, k));

  return result;
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Creates a mapping between moment index n and the index pair (k,l).
 */
vector<Moment> Galerkin_Ordinate_Space::compute_n2lk_3D_(Quadrature_Class /*unused*/,
                                                         unsigned const N) {
  vector<Moment> result;

  // Choose: l= 0, ..., N-1, k = -l, ..., l
  for (unsigned ell = 0; ell < N; ++ell)
    for (int k = -1 * static_cast<int>(ell); k <= static_cast<int>(ell); ++k)
      result.emplace_back(Moment(ell, k));

  // Add ell=N and k<0
  {
    unsigned ell(N);
    for (int k(-1 * static_cast<int>(ell)); k < 0; ++k)
      result.emplace_back(Moment(ell, k));
  }

  // Add ell=N, k>0, k odd
  {
    int ell(N);
    for (int k = 1; k <= ell; k += 2)
      result.emplace_back(Moment(ell, k));
  }

  // Add ell=N+1 and k<0, k even
  {
    unsigned ell(N + 1);
    for (int k(-1 * static_cast<int>(ell) + 1); k < 0; k += 2)
      result.emplace_back(Moment(ell, k));
  }

  return result;
}

//------------------------------------------------------------------------------------------------//
/*!
 * \param dimension Dimension of the physical problem space (1, 2, or 3)
 * \param geometry Geometry of the physical problem space (spherical, axisymmetric, Cartesian)
 * \param ordinates Set of ordinate directions
 * \param quadrature_class Class of the quadrature used to generate the ordinate set. At present,
 *             only TRIANGLE_QUADRATURE is supported.
 * \param sn_order Order of the quadrature. This is equal to the number of levels for triangular and
 *             square quadratures.
 * \param expansion_order Expansion order of the desired scattering moment space.
 * \param method Enum value specifying the desired quadrature type.
 * \param extra_starting_directions Add extra directions to each level set. In most geometries, an
 *             additional ordinate is added that is opposite in direction to the starting
 *             direction. This is used to implement reflection exactly in curvilinear
 *             coordinates. In 1D spherical, that means an additional angle is added at mu=1. In
 *             axisymmetric, that means additional angles are added that are oriented opposite to
 *             the incoming starting direction on each level.
 * \param ordering Ordering into which to sort the ordinates.
 */
Galerkin_Ordinate_Space::Galerkin_Ordinate_Space(
    unsigned const dimension, rtt_mesh_element::Geometry const geometry,
    vector<Ordinate> const &ordinates, Quadrature_Class quadrature_class, unsigned sn_order,
    unsigned const expansion_order, QIM const method, bool const extra_starting_directions,
    Ordering const ordering)
    : Ordinate_Space(dimension, geometry, ordinates, expansion_order, extra_starting_directions,
                     ordering),
      method_(method), D_(), M_() {
  Require(dimension > 0 && dimension < 4);
  Require(geometry != rtt_mesh_element::Geometry::END_GEOMETRY);
  Require(sn_order > 0 && sn_order % 2 == 0);
  Require(method_ == GQ1 || method_ == GQ2 || method_ == GQF);

  // May be relaxed in the future
  Require(expansion_order <= sn_order);
  Require(quadrature_class == TRIANGLE_QUADRATURE || dimension == 1);

  // Creates the moment index map
  compute_moments_(quadrature_class, sn_order);

  // Creates the moment-to-discrete and discrete-to-moment operators
  compute_operators();

  Ensure(check_class_invariants());
}

//------------------------------------------------------------------------------------------------//
bool Galerkin_Ordinate_Space::check_class_invariants() const {
  return D_.size() == ordinates().size() * this->moments().size() &&
         M_.size() == ordinates().size() * this->moments().size();
}

//------------------------------------------------------------------------------------------------//
QIM Galerkin_Ordinate_Space::quadrature_interpolation_model() const {
  Check(method_ == GQ1 || method_ == GQ2 || method_ == GQF);
  return method_;
}

//------------------------------------------------------------------------------------------------//
/*!
 * In the future, this function will allow the client to specify the maximum order to include, but
 * for now, we include all full orders, leaving out any Galerkin augments.
 */
vector<double> Galerkin_Ordinate_Space::D() const {
  size_t const number_of_ordinates = ordinates().size();
  unsigned const number_of_moments = this->number_of_moments();

  vector<double> Result(number_of_ordinates * number_of_moments);

  for (unsigned a = 0; a < number_of_ordinates; ++a) {
    for (unsigned m = 0; m < number_of_moments; ++m) {
      Result[a + number_of_ordinates * m] = D_[a + number_of_ordinates * m];
    }
  }
  return Result;
}
//------------------------------------------------------------------------------------------------//
/*!
 * In the future, this function will allow the client to specify the maximum order to include, but
 * for now, we include all full orders, leaving out any Galerkin augments.
 */
vector<double> Galerkin_Ordinate_Space::M() const {
  size_t const number_of_ordinates = ordinates().size();
  unsigned const number_of_moments = this->number_of_moments();
  size_t const total_number_of_moments = this->moments().size();

  vector<double> Result(number_of_ordinates * number_of_moments);

  for (unsigned a = 0; a < number_of_ordinates; ++a) {
    for (unsigned m = 0; m < number_of_moments; ++m) {
      Result[m + number_of_moments * a] = M_[m + total_number_of_moments * a];
    }
  }
  return Result;
}

//------------------------------------------------------------------------------------------------//
void Galerkin_Ordinate_Space::compute_operators() {

  rtt_mesh_element::Geometry const geometry(this->geometry());

  vector<Ordinate> &ordinates(this->ordinates());
  size_t const numOrdinates(ordinates.size());

  vector<Ordinate> cartesian_ordinates;
  vector<unsigned> indexes;

  // fill cartesian_ordinates
  unsigned count(0);
  for (unsigned i = 0; i < numOrdinates; ++i) {
    if (std::abs(ordinates[i].wt()) > std::numeric_limits<decltype(ordinates[i].wt())>::min()) {
      cartesian_ordinates.emplace_back(ordinates[i]);
      indexes.emplace_back(count++);
    } else
      indexes.emplace_back(0);
  }

  size_t const numCartesianOrdinates(cartesian_ordinates.size());
  size_t const numMoments(this->moments().size());

  // create Cartesian SN operators
  vector<double> cartesian_M_SN(compute_M_SN(cartesian_ordinates));

  vector<double> cartesian_M;
  vector<double> cartesian_D;

  if (method_ == 1 || method_ == 3) {
    // ------------------------------------------------------------------------
    // invert the (m x n) moment-to-discrete matrix M to compute the discrete-to-moment matrix D
    // ------------------------------------------------------------------------

    cartesian_M.swap(cartesian_M_SN);
    Check(numMoments < UINT_MAX);
    Check(numCartesianOrdinates < UINT_MAX);
    cartesian_D = compute_inverse(static_cast<unsigned>(numMoments),
                                  static_cast<unsigned>(numCartesianOrdinates), cartesian_M);

    // set cartesian_ordinate weights to the first row of D

    for (unsigned i = 0; i < numCartesianOrdinates; ++i) {
      cartesian_ordinates[i].set_wt(cartesian_D[i + 0 * numCartesianOrdinates]);
    }

    // and reset ordinate weights to the first row of D

    vector<Ordinate> &lordinates(this->ordinates());
    for (unsigned i = 0; i < numOrdinates; ++i) {
      if (std::abs(lordinates[i].wt()) > std::numeric_limits<decltype(lordinates[i].wt())>::min()) {
        lordinates[i].set_wt(cartesian_D[indexes[i] + 0 * numCartesianOrdinates]);
      }
    }
  } else if (method_ == 2) {
    // first get new ordinate weights from the usual GQ method, needed to accurately integrate all
    // moments

    // This branch is not tested
    Insist(false, "This branch commented out because there are no supporting "
                  "unit tests");
  } else
    Insist(false, "Could not identify Galerkin Quadrature method.");

  // store the final form of the operators in M_ and D_

  if (geometry == rtt_mesh_element::Geometry::CARTESIAN) {
    M_ = cartesian_M;
    D_ = cartesian_D;
  } else // augment the Cartesian operators for the zero-weight starting directions then store
  {
    M_ = augment_M(indexes, cartesian_M);

    Check(numCartesianOrdinates < UINT_MAX);
    D_ = augment_D(indexes, static_cast<unsigned>(numCartesianOrdinates), cartesian_D);
  }

  for (unsigned n = 0; n < numMoments; ++n) {
    unsigned const ell(moments()[n].L());
    int const k(moments()[n].M());

    std::cout << " moment " << n << "     l = " << ell << " k = " << k << std::endl;
  }
}

//------------------------------------------------------------------------------------------------//
// Augment the matrix for curvilinear coordinates
vector<double> Galerkin_Ordinate_Space::augment_D(vector<unsigned> const &indexes,
                                                  unsigned const numCartesianOrdinates,
                                                  vector<double> const &D) {
  vector<Ordinate> const &ordinates(this->ordinates());
  size_t const numOrdinates(ordinates.size());

  size_t const numMoments(this->moments().size());

  Check(indexes.size() == numOrdinates);

  vector<double> D_new(numMoments * numOrdinates, 0);

  for (unsigned m = 0; m < numOrdinates; ++m) {
    for (unsigned n = 0; n < numMoments; ++n) {
      if (std::abs(ordinates[m].wt()) > std::numeric_limits<decltype(ordinates[m].wt())>::min()) {
        D_new[m + n * numOrdinates] = D[indexes[m] + n * numCartesianOrdinates];
      }
    }
  }
  return D_new;
}

//------------------------------------------------------------------------------------------------//
// Augment the matrix for curvilinear coordinates
vector<double> Galerkin_Ordinate_Space::augment_M(vector<unsigned> const &indexes,
                                                  vector<double> const &M) {
  using rtt_sf::Ylm;

  vector<Ordinate> const &ordinates(this->ordinates());
  size_t const numOrdinates(ordinates.size());

  vector<Moment> const &n2lk(this->moments());
  size_t const numMoments(n2lk.size());

  double const sumwt(norm());

  Check(indexes.size() == numOrdinates);

  vector<double> M_new(numMoments * ordinates.size(), 0);

  for (unsigned n = 0; n < numMoments; ++n) {
    unsigned const ell(n2lk[n].L());
    int const k(n2lk[n].M());

    for (unsigned m = 0; m < numOrdinates; ++m) {
      if (std::abs(ordinates[m].wt()) > std::numeric_limits<decltype(ordinates[m].wt())>::min()) {
        M_new[n + m * numMoments] = M[n + indexes[m] * numMoments];
      } else {
        double mu(ordinates[m].mu());
        double eta(ordinates[m].eta());
        double xi(ordinates[m].xi());

        double phi(compute_azimuthalAngle(mu, xi));
        M_new[n + m * numMoments] = Ylm(ell, k, eta, phi, sumwt);
      }
    }
  }
  return M_new;
}

//------------------------------------------------------------------------------------------------//
vector<double> Galerkin_Ordinate_Space::compute_M_SN(vector<Ordinate> const &ordinates) {
  using rtt_sf::Ylm;

  rtt_mesh_element::Geometry const geometry(this->geometry());
  unsigned const dim(dimension());

  vector<Moment> const &n2lk(moments());
  size_t const numMoments(n2lk.size());
  size_t const numOrdinates(ordinates.size());
  double const sumwt(norm());

  // resize the M matrix.
  std::vector<double> M(numMoments * numOrdinates);

  //    double polar, azimuthal;
  for (unsigned n = 0; n < numMoments; ++n) {
    unsigned const ell(n2lk[n].L());
    int const k(n2lk[n].M());

    for (unsigned m = 0; m < numOrdinates; ++m) {
      if (dim == 1 &&
          geometry != rtt_mesh_element::Geometry::AXISYMMETRIC) // 1D mesh, 1D quadrature
      {
        double mu(ordinates[m].mu());
        M[n + m * numMoments] = Ylm(ell, k, mu, 0.0, sumwt);
      } else {
        double mu(ordinates[m].mu());
        double eta(ordinates[m].eta());
        double xi(ordinates[m].xi());

        if (geometry == rtt_mesh_element::Geometry::AXISYMMETRIC) {
          // R-Z coordinate system
          //
          // It is important to remember here that the positive mu axis points to the left and the
          // positive eta axis points up, when the unit sphere is projected on the plane of the mu-
          // and eta-axis in R-Z. In this case, phi is measured from the mu-axis counterclockwise.
          //
          // This accounts for the fact that the azimuthal angle is discretized on levels of the
          // xi-axis, making the computation of the azimuthal angle here consistent with the
          // discretization by using the eta and mu ordinates to define phi.

          double phi(compute_azimuthalAngle(mu, xi));
          M[n + m * numMoments] = Ylm(ell, k, eta, phi, sumwt);
        } else if (geometry == rtt_mesh_element::Geometry::CARTESIAN) {
          // X-Y coordinate system
          //
          // In order to make the harmonic trial space is correctly oriented with respect to the
          // moments chosen, the value of xi and eta are swapped.

          double phi(compute_azimuthalAngle(mu, eta));
          M[n + m * numMoments] = Ylm(ell, k, xi, phi, sumwt);
        }
      }
    } // ordinate loop
  }   // moment loop

  return M;
}

//------------------------------------------------------------------------------------------------//
vector<double> Galerkin_Ordinate_Space::compute_inverse(unsigned const m, unsigned const n,
                                                        vector<double> const &Ain) {
  // Invert an (m x n) matrix A

  Insist(!Ain.empty(), "The GQ ordinate space computation for D requires "
                       "that M be available.");

  Insist(n == m, "Matrix must be square.");

  std::vector<double> A(Ain);
  std::vector<double> B(m * n);

  // Create GSL matrix views of A and B.  LU will get a copy of M.  This matrix will be decomposed
  // into LU.
  gsl_matrix_view gsl_A = gsl_matrix_view_array(&A[0], m, n);
  gsl_matrix_view gsl_B = gsl_matrix_view_array(&B[0], n, m);

  // Create some local space for the permutation matrix.
  gsl_permutation *p = gsl_permutation_alloc(m);

  // Store information about sign changes in this variable.
  int signum(0);

  // Factorize the square matrix M into the LU decomposition PM = LU.  On output the diagonal and
  // upper triangular part of the input matrix M contain the matrix U.  The lower triangular part of
  // the input matrix (excluding the diagonal) contains L. The diagonal elements of L are unity, and
  // are not stored.
  //
  // The permutation matrix P is encoded in the permutation p.  The j-th column of the matrix P is
  // given by the k-th column of the identity, where k=p[j] the j-th element of the permutation
  // vector.  The sign of the permutation is given by signum.  It has the value \f$ (-1)^n \f$,
  // where n is the number of interchanges in the permutation.
  //
  // The algorithm used in the decomposition is Gaussian Elimination with partial pivoting (Golub &
  // Van Loan, Matrix Computations, Algorithm 3.4.1).

  // Store the LU decomposition in the matrix A.
  Remember(int result =) gsl_linalg_LU_decomp(&gsl_A.matrix, p, &signum);
  Check(result == 0);

  // Compute the inverse of the matrix LU from its LU decomposition (LU,p), storing the results in
  // the matrix B.  The inverse is computed by solving the system (LU) x = b for each column of the
  // identity matrix.

  Remember(result =) gsl_linalg_LU_invert(&gsl_A.matrix, p, &gsl_B.matrix);

  Check(result == 0);

  // Free the space reserved for the permutation matrix.
  gsl_permutation_free(p);

  return B;
}

} // end namespace rtt_quadrature

//------------------------------------------------------------------------------------------------//
// end of Galerkin_Ordinate_Space.cc
//------------------------------------------------------------------------------------------------//
