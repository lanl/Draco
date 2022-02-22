//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   quadrature/Ordinate_Set_Factory.cc
 * \author Allan Wollaber
 * \date   Mon Mar  7 10:42:56 EST 2016
 * \brief  Implementation file for the class rtt_quadrature::Ordinate_Set_Factory.
 * \note   Copyright (C) 2016-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "Ordinate_Set_Factory.hh"
#include "Gauss_Legendre.hh"
#include "Level_Symmetric.hh"
#include "Lobatto.hh"
#include "Product_Chebyshev_Legendre.hh"
#include "Quadrature.hh"
#include "Quadrature_Interface.hh"
#include "Square_Chebyshev_Legendre.hh"
#include "Tri_Chebyshev_Legendre.hh"
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <numeric>

namespace rtt_quadrature {

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Primary service method to create an Ordinate_Set
 *
 * \return A smart pointer to the correctly instantiated Ordinate_Set
 *
 *  All of the "magic numbers" in the quadrature_data struct get interpreted here to build the
 *  Ordinate_Set object.
 *
 *  This class/method combination would be better suited as a factory method (and not a factory
 *  class) within Ordinate_Set, but this would introduce a cyclic dependency between Quadrature and
 *  Ordinate_Set. Alternatively, it could be a factory method in Quadrature, but we'd then need an
 *  additional factory method to create the Quadrature object itself. As seen below, the derived
 *  Quadrature object need only exist temporarily in order to call "create_ordinate_set".  This
 *  solution avoids both of those pitfalls, but it is not ideal and could be refactored into one or
 *  the other classes if their design ever changes.
 */
std::shared_ptr<Ordinate_Set> Ordinate_Set_Factory::get_Ordinate_Set() const {

  using rtt_mesh_element::Geometry;

  bool add_starting_directions = false;
  bool add_extra_directions = false;

  rtt_mesh_element::Geometry geometry;

  // Find the geometry
  switch (quad_.geometry) {
  case 0:
    geometry = rtt_mesh_element::Geometry::CARTESIAN;
    break;

  case 1:
    geometry = rtt_mesh_element::Geometry::AXISYMMETRIC;
    add_starting_directions = true;
    break;

  case 2:
    geometry = rtt_mesh_element::Geometry::SPHERICAL;
    add_starting_directions = true;
    break;

  default:
    Insist(false, "Unrecognized Geometry");
    geometry = rtt_mesh_element::Geometry::CARTESIAN;
  }

  std::shared_ptr<Ordinate_Set> ordinate_set;

  if (quad_.dimension == 1) { // 1D quadratures

    if (quad_.type == 0) {
      Gauss_Legendre quadrature(quad_.order);
      ordinate_set = quadrature.create_ordinate_set(1, geometry,
                                                    1.0, // norm,
                                                    add_starting_directions, add_extra_directions,
                                                    Ordinate_Set::LEVEL_ORDERED);
    } else if (quad_.type == 1) {
      Lobatto quadrature(quad_.order);
      ordinate_set = quadrature.create_ordinate_set(1, geometry,
                                                    1.0, // norm,
                                                    add_starting_directions, add_extra_directions,
                                                    Ordinate_Set::LEVEL_ORDERED);
    }
  } else if (quad_.dimension == 2) { // 2D quadratures
    if (quad_.type == 0) {
      Level_Symmetric quadrature(quad_.order);
      ordinate_set = quadrature.create_ordinate_set(2, geometry,
                                                    1.0, // norm,
                                                    add_starting_directions, add_extra_directions,
                                                    Ordinate_Set::LEVEL_ORDERED);
    } else if (quad_.type == 1) {
      Tri_Chebyshev_Legendre quadrature(quad_.order);
      ordinate_set = quadrature.create_ordinate_set(2, geometry,
                                                    1.0, // norm,
                                                    add_starting_directions, add_extra_directions,
                                                    Ordinate_Set::LEVEL_ORDERED);
    } else if (quad_.type == 2) {
      Square_Chebyshev_Legendre quadrature(quad_.order);
      ordinate_set = quadrature.create_ordinate_set(2, geometry,
                                                    1.0, // norm,
                                                    add_starting_directions, add_extra_directions,
                                                    Ordinate_Set::LEVEL_ORDERED);
    } else if (quad_.type == 3) {
      Product_Chebyshev_Legendre quadrature(quad_.order, quad_.azimuthal_order);
      ordinate_set = quadrature.create_ordinate_set(2, geometry,
                                                    1.0, // norm,
                                                    add_starting_directions, add_extra_directions,
                                                    Ordinate_Set::LEVEL_ORDERED);
    }
  } else {
    Insist(quad_.dimension > 0 && quad_.dimension < 3, "Expected 1 <= quad_.dimension < 3.");
  }

  Ensure(ordinate_set);
  return ordinate_set;
}

} // end namespace rtt_quadrature

//------------------------------------------------------------------------------------------------//
// end of quadrature/Ordinate_Set_Factory.cc
//------------------------------------------------------------------------------------------------//
