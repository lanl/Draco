//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   quadrature/Interval_Quadrature.cc
 * \author Kelly Thompson
 * \date   Tue Feb 22 10:21:50 2000
 * \brief  A class representing an interval Gauss-Legendre quadrature set.
 * \note   Copyright (C) 2012-2022 Triad National Security, LLC., All rights reserved.
 *         reserved. */
//------------------------------------------------------------------------------------------------//

#include "Interval_Quadrature.hh"

namespace rtt_quadrature {
using namespace std;

Interval_Quadrature::Interval_Quadrature(unsigned const sn_order) : Quadrature(sn_order) {}

//------------------------------------------------------------------------------------------------//
/* virtual */
Quadrature_Class Interval_Quadrature::quadrature_class() const { return INTERVAL_QUADRATURE; }

//------------------------------------------------------------------------------------------------//

bool Interval_Quadrature::has_axis_assignments() const {
  return false; // cannot override default assignments
}

//------------------------------------------------------------------------------------------------//
/* virtual */
vector<Ordinate> Interval_Quadrature::create_ordinates_(unsigned const /*dimension*/,
                                                        Geometry const geometry, double const norm,
                                                        unsigned const /*mu_axis*/,
                                                        unsigned /*eta_axis*/,
                                                        bool const include_starting_directions,
                                                        bool const include_extra_directions) const {
  vector<Ordinate> Result = create_level_ordinates_(norm);

  // add any starting or extra directions

  add_1D_starting_directions_(geometry, include_starting_directions, include_extra_directions,
                              Result);

  return Result;
}

//------------------------------------------------------------------------------------------------//
/* virtual */
vector<Ordinate> Interval_Quadrature::create_ordinates_(unsigned const /*dimension*/,
                                                        Geometry const geometry, double const norm,
                                                        bool const include_starting_directions,
                                                        bool const include_extra_directions) const {
  return Interval_Quadrature::create_ordinates_(1, // can only be 1-D
                                                geometry, norm,
                                                0, // can only be aligned with mu
                                                0, include_starting_directions,
                                                include_extra_directions);
}

} // end namespace rtt_quadrature

//------------------------------------------------------------------------------------------------//
// end of quadrature/Interval_Quadrature.cc
//------------------------------------------------------------------------------------------------//
