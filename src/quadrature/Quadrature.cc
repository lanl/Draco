//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   quadrature/Quadrature.cc
 * \author Kelly Thompson
 * \date   Tue Feb 22 15:38:56 2000
 * \brief  Copyright (C) 2010-2021 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "Quadrature.hh"
#include "Galerkin_Ordinate_Space.hh"
#include "Sn_Ordinate_Space.hh"
#include <algorithm>

namespace rtt_quadrature {

using std::make_shared;
using Ordering = Ordinate_Set::Ordering;
using rtt_dsxx::soft_equiv;

//------------------------------------------------------------------------------------------------//
/*!
 * Create a set of ordinates from the Quadrature.
 *
 * \param dimension Dimension of the problem.
 * \param geometry Geometry of the problem.
 * \param norm Norm for the ordinate weights.
 * \param mu_axis Which spatial axis maps to the mu direction cosine?
 * \param eta_axis Which spatial axis maps to the eta direction cosine?
 * \param include_starting_directions Should starting directions be included in the ordinate set for
 *          each level? This argument is ignored for Cartesian geometries.
 * \param include_extra_directions Should extra starting directions be included in the ordinate set
 *          for each level?
 */
std::vector<Ordinate> Quadrature::create_ordinates(unsigned const dimension,
                                                   Geometry const geometry, double const norm,
                                                   unsigned const mu_axis, unsigned const eta_axis,
                                                   bool const include_starting_directions,
                                                   bool const include_extra_directions) const {
  Require(dimension > 0 && dimension < 4);
  Require(norm > 0.0);
  Require(mu_axis < 3 && eta_axis < 3);
  Require(dimension > 1 || mu_axis == 0);
  Require(dimension == 1 || mu_axis != eta_axis);
  Require(dimension == 1 || quadrature_class() != INTERVAL_QUADRATURE);

  std::vector<Ordinate> Result =
      create_ordinates_(dimension, geometry, norm, mu_axis, eta_axis, include_starting_directions,
                        include_extra_directions);

  return Result;
}

//------------------------------------------------------------------------------------------------//
/*!
 * Create a set of ordinates from the Quadrature.
 *
 * \param dimension Dimension of the problem.
 * \param geometry Geometry of the problem.
 * \param norm Norm for the ordinate weights.
 * \param include_starting_directions Should starting directions be included in the ordinate set for
 *          each level? This argument is ignored for Cartesian geometries.
 * \param include_extra_directions Should extra starting directions be included in the ordinate set
 *          for each level?
 */
std::vector<Ordinate> Quadrature::create_ordinates(unsigned const dimension,
                                                   Geometry const geometry, double const norm,
                                                   bool const include_starting_directions,
                                                   bool const include_extra_directions) const {
  Require(dimension > 0 && dimension < 4);
  Require(norm > 0.0);
  Require(dimension == 1 || quadrature_class() != INTERVAL_QUADRATURE);

  std::vector<Ordinate> Result = create_ordinates_(
      dimension, geometry, norm, include_starting_directions, include_extra_directions);

  return Result;
}

//------------------------------------------------------------------------------------------------//
/*!
 * Create an ordinate set from the Quadrature.
 *
 * \param dimension Dimension of the problem.
 * \param geometry Geometry of the problem.
 * \param norm Norm for the ordinate weights.
 * \param include_starting_directions Should starting directions be included in the ordinate set for
 *          each level? This argument is ignored for Cartesian geometries.
 * \param include_extra_directions Should extra starting directions be included in the ordinate set
 *          for each level?
 * \param ordering What ordering should be imposed on the ordinates?
 */
std::shared_ptr<Ordinate_Set>
Quadrature::create_ordinate_set(unsigned const dimension, Geometry const geometry,
                                double const norm, bool const include_starting_directions,
                                bool const include_extra_directions,
                                Ordering const ordering) const {
  Require(dimension > 0 && dimension < 4);
  Require(norm > 0.0);
  Require(dimension == 1 || quadrature_class() != INTERVAL_QUADRATURE);

  std::vector<Ordinate> ordinates = create_ordinates_(
      dimension, geometry, norm, include_starting_directions, include_extra_directions);

  std::shared_ptr<Ordinate_Set> Result(new Ordinate_Set(dimension, geometry, ordinates,
                                                        include_starting_directions,
                                                        include_extra_directions, ordering));

  return Result;
}

//------------------------------------------------------------------------------------------------//
/* protected */
void Quadrature::add_1D_starting_directions_(Geometry const geometry,
                                             bool const add_starting_directions,
                                             bool const add_extra_starting_directions,
                                             std::vector<Ordinate> &ordinates) const {
  // Add starting directions if necessary

  if (add_starting_directions) {
    if (geometry == rtt_mesh_element::Geometry::SPHERICAL) {
      // insert mu=-1 starting direction
      auto a = ordinates.begin();
      a = ordinates.insert(a, Ordinate(-1.0, 0.0, 0.0, 0.0));

      // insert mu=1 starting direction
      if (add_extra_starting_directions)
        ordinates.emplace_back(Ordinate(1.0, 0.0, 0.0, 0.0));
    } else if (geometry == rtt_mesh_element::Geometry::AXISYMMETRIC) {
      Insist(false, "should not be reached");
    }
  }
}

//---------------------------------------------------------------------------------------//
/* protected */
void Quadrature::add_2D_starting_directions_(Geometry const geometry,
                                             bool const add_starting_directions,
                                             bool const add_extra_starting_directions,
                                             std::vector<Ordinate> &ordinates) const {
  // Add starting directions if necessary

  if (add_starting_directions) {
    if (geometry == rtt_mesh_element::Geometry::AXISYMMETRIC) {
      std::sort(ordinates.begin(), ordinates.end(), Ordinate_Set::level_compare);

      // Define an impossible value for a direction cosine.  We use this to simplify the logic of
      // determining when we are at the head of a new level set.

      double const SENTINEL_COSINE = 2.0;

      // Insert the supplemental ordinates.

      double eta = -SENTINEL_COSINE;

      for (auto a = ordinates.begin(); a != ordinates.end(); ++a) {
        double const old_eta = eta;
        eta = a->eta();
        if (!soft_equiv(eta, old_eta))
        // We are at the start of a new level.  Insert the starting ordinate.  This has xi==0 and mu
        // determined by the normalization condition.
        {
          Check(1.0 - eta * eta >= 0.0);

          // insert mu < 0
          a = ordinates.insert(a, Ordinate(-sqrt(1.0 - eta * eta), eta, 0.0, 0.0));

          // insert mu > 0
          if (add_extra_starting_directions)
            if (a != ordinates.begin())
              a = ordinates.insert(a, Ordinate(sqrt(1.0 - old_eta * old_eta), old_eta, 0.0, 0.0));
        }
      }

      // insert mu > 0 on the final level
      if (add_extra_starting_directions)
        ordinates.emplace_back(Ordinate(sqrt(1.0 - eta * eta), eta, 0.0, 0.0));
    }
  }
}

//------------------------------------------------------------------------------------------------//
void Quadrature::map_axes_(unsigned const mu_axis, unsigned const eta_axis, std::vector<double> &mu,
                           std::vector<double> &eta, std::vector<double> &xi) const {
  if (mu_axis == 0) {
    if (eta_axis == 1) {
      // no action needed -- defaults in place
    } else {
      eta.swap(xi);
    }
  } else if (mu_axis == 1) {
    if (eta_axis == 0) {
      mu.swap(eta);
    } else {
      eta.swap(xi);
      mu.swap(xi);
    }
  } else {
    if (eta_axis == 0) {
      mu.swap(eta);
      mu.swap(xi);
    } else {
      mu.swap(xi);
    }
  }
}

//------------------------------------------------------------------------------------------------//
/*!
 * Create an Ordinate_Space from the Quadrature.
 *
 * \param dimension Dimension of the problem.
 * \param geometry Geometry of the problem.
 * \param moment_expansion_order Expansion order in moment space
 * \param mu_axis Which spatial axis maps to the mu direction cosine?
 * \param eta_axis Which spatial axis maps to the eta direction cosine?
 * \param include_extra_directions Should extra starting directions be included in the ordinate set
 *          for each level?
 * \param ordering What ordering should be imposed on the ordinates?
 * \param qim What interpolation model should be used to generate the moment space?
 */
std::shared_ptr<Ordinate_Space>
Quadrature::create_ordinate_space(unsigned const dimension, Geometry const geometry,
                                  unsigned const moment_expansion_order, unsigned const mu_axis,
                                  unsigned const eta_axis, bool const include_extra_directions,
                                  Ordering const ordering, QIM const qim) const {
  Require(dimension > 0 && dimension < 4);
  Require(mu_axis < 3 && eta_axis < 3);
  Require(dimension > 1 || mu_axis == 0);
  Require(dimension == 1 || mu_axis != eta_axis);
  Require(dimension == 1 || quadrature_class() != INTERVAL_QUADRATURE);

  std::vector<Ordinate> ordinates = create_ordinates(dimension, geometry,
                                                     1.0, // hardwired norm
                                                     mu_axis, eta_axis,
                                                     true, // include starting directions
                                                     include_extra_directions);

  std::shared_ptr<Ordinate_Space> Result;
  switch (qim) {
  case SN:
    Result = std::make_shared<Sn_Ordinate_Space>(
        dimension, geometry, ordinates, moment_expansion_order, include_extra_directions, ordering);
    break;

  case GQ1:
    Result = std::make_shared<Galerkin_Ordinate_Space>(
        dimension, geometry, ordinates, quadrature_class(), number_of_levels(),
        moment_expansion_order, GQ1, include_extra_directions, ordering);
    break;

  case GQ2:
    Result = std::make_shared<Galerkin_Ordinate_Space>(
        dimension, geometry, ordinates, quadrature_class(), number_of_levels(),
        moment_expansion_order, GQ2, include_extra_directions, ordering);
    break;

  default:
    Insist(false, "bad case");
  }

  return Result;
}

//------------------------------------------------------------------------------------------------//
/*!
 * Create an angle operator from the Quadrature.
 *
 * \param dimension Dimension of the problem.
 * \param geometry Geometry of the problem.
 * \param moment_expansion_order Expansion order in moment space. If negative, the moment space is
 *          not needed.
 * \param include_extra_directions Should extra starting directions be included in the ordinate set
 *          for each level?
 * \param ordering What ordering should be imposed on the ordinates?
 * \param qim What interpolation model should be used to generate the moment space?
 */
std::shared_ptr<Ordinate_Space> Quadrature::create_ordinate_space(
    unsigned const dimension, Geometry const geometry, int const moment_expansion_order,
    bool const include_extra_directions, Ordering const ordering, QIM const qim) const {
  Require(dimension > 0 && dimension < 4);
  Require(dimension == 1 || quadrature_class() != INTERVAL_QUADRATURE);
  Require(qim == SN || qim == GQ1 || qim == GQ2 || qim == GQF);
  Require(qim == SN || moment_expansion_order >= 0);

  std::vector<Ordinate> ordinates =
      create_ordinates(dimension, geometry,
                       1.0, // hardwired norm
                       geometry != rtt_mesh_element::Geometry::CARTESIAN,
                       // include starting directions if curvilinear
                       include_extra_directions);

  std::shared_ptr<Ordinate_Space> Result;
  if (qim == SN)
    Result = make_shared<Sn_Ordinate_Space>(dimension, geometry, ordinates, moment_expansion_order,
                                            include_extra_directions, ordering);
  else
    Result = make_shared<Galerkin_Ordinate_Space>(
        dimension, geometry, ordinates, quadrature_class(), number_of_levels(),
        moment_expansion_order, qim, include_extra_directions, ordering);

  return Result;
}

//------------------------------------------------------------------------------------------------//
bool Quadrature::is_open_interval() const {
  // The great majority are. Lobatto and certain cases of General Octant are at present our only
  // exceptions.
  return true;
}

} // end namespace rtt_quadrature

//------------------------------------------------------------------------------------------------//
// end of Quadrature.cc
//------------------------------------------------------------------------------------------------//
