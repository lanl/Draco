//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   quadrature/Sn_Ordinate_Space.hh
 * \author Kent Budge
 * \date   Mon Mar 26 16:11:19 2007
 * \brief  Definition of class Sn_Ordinate_Space
 * \note   Copyright (C) 2012-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef quadrature_Sn_Ordinate_Space_hh
#define quadrature_Sn_Ordinate_Space_hh

#include "Ordinate_Space.hh"

namespace rtt_quadrature {

//================================================================================================//
/*!
 * \class Sn_Ordinate_Space
 *
 * \brief Represents ordinate operators for a conventional Sn moment space.
 *
 * The moment space contains all moments up to the specified scattering order, and the moment to
 * discrete and discrete to moment operators are calculated in a straightforward manner from the Ylm
 * and the weight associated with each ordinate direction.
 */
//================================================================================================//

class Sn_Ordinate_Space : public Ordinate_Space {
public:
  // CREATORS

  //! Specify the ordinate quadrature with defaults.
  Sn_Ordinate_Space(unsigned dimension, rtt_mesh_element::Geometry geometry,
                    std::vector<Ordinate> const &ordinates, int expansion_order,
                    bool extra_starting_directions = false, Ordering ordering = LEVEL_ORDERED);

  // ACCESSORS

  bool check_class_invariants() const;

  // SERVICES

  QIM quadrature_interpolation_model() const override;

  //! Return the discrete to moment transform matrix
  std::vector<double> D() const override;

  //! Return the moment to discrete transform matrix
  std::vector<double> M() const override;

protected:
  // IMPLEMENTATION

  std::vector<Moment> compute_n2lk_1D_(Quadrature_Class /*unused*/, unsigned sn_order) override;
  std::vector<Moment> compute_n2lk_1Da_(Quadrature_Class /*unused*/, unsigned sn_order) override;
  std::vector<Moment> compute_n2lk_2D_(Quadrature_Class /*unused*/, unsigned sn_order) override;
  std::vector<Moment> compute_n2lk_2Da_(Quadrature_Class /*unused*/, unsigned sn_order) override;
  std::vector<Moment> compute_n2lk_3D_(Quadrature_Class /*unused*/, unsigned sn_order) override;

private:
  // IMPLEMENTATION

  void compute_M();
  void compute_D();

  // DATA

  std::vector<double> D_; //!< Discrete to Moment matrix
  std::vector<double> M_; //!< Moment to Discrete matrix
};

} // end namespace rtt_quadrature

#endif // quadrature_Sn_Ordinate_Space_hh

//------------------------------------------------------------------------------------------------//
// end of quadrature/Sn_Ordinate_Space.hh
//------------------------------------------------------------------------------------------------//
