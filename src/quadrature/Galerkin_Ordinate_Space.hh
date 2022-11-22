//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   quadrature/Galerkin_Ordinate_Space.hh
 * \author Kent Budge
 * \date   Mon Mar 26 16:11:19 2007
 * \brief  Definition of class Galerkin_Ordinate_Space
 * \note   Copyright (C) 2012-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef quadrature_Galerkin_Ordinate_Space_hh
#define quadrature_Galerkin_Ordinate_Space_hh

#include "Ordinate_Space.hh"

namespace rtt_quadrature {
using std::ostream;

//================================================================================================//
/*!
 * \class Galerkin_Ordinate_Space
 * \brief Represents ordinate operators for a Galerkin moment space.
 *
 * The moment space contains all moments (that are not identically zero due to symmetry) up to the
 * specified scattering order, but the moment to discrete operator M and discrete to moment operator
 * D are computed as if enough additional higher moments are included in the moment space to make D
 * and M square. The higher moment terms are then discarded, but the non-square D and M retain the
 * property that DM is the identity. This stabilizes the moment to discrete and discrete to moment
 * operations at high scattering orders.
 *
 * When the additional moments are added, the SN quadrature order is provided, and additional
 * moments added based on the assumption that we are using triangular quadrature sets. If an
 * expansion order L < N-1 is requested, both D and M will be appropriately truncated. If an
 * expansion order L > N-1 is requested, the computation cannot proceed and an exception thrown.
 *
 * \todo For triangular sets, we assume L = N-1. The other two cases, L > N-1 and L < N-1 need to be
 *       handled.
 *
 * \todo The Galerkin quadrature is only currently implemented for triangular quadratures, such that
 *       it will be necessary to extend the method to include square (product) quadratures.
 */
//================================================================================================//

class Galerkin_Ordinate_Space : public Ordinate_Space {
public:
  // CREATORS

  //! Specify the ordinate quadrature with defaults.
  Galerkin_Ordinate_Space(unsigned const dimension, rtt_mesh_element::Geometry const geometry,
                          std::vector<Ordinate> const &ordinates, Quadrature_Class quadrature_class,
                          unsigned sn_order, unsigned const expansion_order, QIM const method,
                          bool const extra_starting_directions = false,
                          Ordering const ordering = LEVEL_ORDERED);

  // ACCESSORS

  bool check_class_invariants() const;

  // SERVICES

  QIM quadrature_interpolation_model() const override;

  //! Return the discrete to moment transform matrix
  std::vector<double> D() const override;

  //! Return the moment to discrete transform matrix
  std::vector<double> M() const override;

  /*! Prune any moments beyond the user-specified expansion order. Such moments are included in
   *  Galerkin methods for purposes of computing the M and D matrices, but are then removed from the
   *  moment space unless the GQF interpolation model has been specified. */
  bool prune() const override { return method_ != GQF; }

  // STATICS

private:
  // IMPLEMENTATION

  std::vector<Moment> compute_n2lk_1D_(Quadrature_Class /*unused*/, unsigned sn_order) override;
  std::vector<Moment> compute_n2lk_1Da_(Quadrature_Class /*unused*/, unsigned sn_order) override;
  std::vector<Moment> compute_n2lk_2D_(Quadrature_Class /*unused*/, unsigned sn_order) override;
  std::vector<Moment> compute_n2lk_2Da_(Quadrature_Class /*unused*/, unsigned sn_order) override;
  std::vector<Moment> compute_n2lk_3D_(Quadrature_Class /*unused*/, unsigned sn_order) override;

private:
  // IMPLEMENTATION

  void compute_operators();

  std::vector<double> compute_M_SN(std::vector<Ordinate> const &ordinates);
  std::vector<double> compute_D_SN(std::vector<Ordinate> const &ordinates,
                                   std::vector<double> const &Min);

  std::vector<double> compute_inverse(unsigned const m, unsigned const n,
                                      std::vector<double> const &Ain);

  std::vector<double> augment_D(std::vector<unsigned> const &indexes,
                                unsigned const numCartesianOrdinates, std::vector<double> const &D);

  std::vector<double> augment_M(std::vector<unsigned> const &indexes, std::vector<double> const &M);

  // DATA
  QIM const method_;
  std::vector<double> D_; //!< Discrete to moment matrix
  std::vector<double> M_; //!< Moment to discrete matrix
};

} // end namespace rtt_quadrature

#endif // quadrature_Galerkin_Ordinate_Space_hh

//------------------------------------------------------------------------------------------------//
// end of quadrature/Galerkin_Ordinate_Space.hh
//------------------------------------------------------------------------------------------------//
