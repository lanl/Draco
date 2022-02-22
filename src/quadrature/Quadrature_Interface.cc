//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   quadrature/Quadrature_Interface.cc
 * \author Jae Chang
 * \date   Tue Jan 27 08:51:19 2004
 * \brief  Quadrature interface definitions
 * \note   Copyright (C) 2015-2022 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

#include "Quadrature_Interface.hh"
#include "Ordinate_Set_Factory.hh"
#include <iostream>

using std::vector;

//! An extern "C" interface to default constructor
void init_quadrature(quadrature_data &quad) { quad = quadrature_data(); }

// Function to check basic validity of quadrature_data
void check_quadrature_validity(const quadrature_data &quad) {
  // We only support 1, 2, and 3D problems
  Insist(quad.dimension > 0 && quad.dimension <= 3, "Quadrature dimension must be 1, 2, or 3");

  if (quad.dimension == 1) {
    Insist(0 <= quad.type && quad.type <= 2, "Quadrature type must be 1 or 2 in 1-D");
  } else if (quad.dimension == 2) {
    Insist(0 <= quad.type && quad.type <= 3, "Quadrature type must be in [0,3] in 2-D");
  }

  // Check for valid order
  Insist(quad.order > 0, "Quadrature order must be positive");

  // There are no checks on azimuthal order since it's not required.

  // There are 3 geometry types (0,1,2) supported
  Insist(0 <= quad.geometry && quad.geometry <= 2, "Invalid geometry in quadrature_data");

  // The "mu" and "weights" entries must not be NULL; others can be
  Insist(quad.mu != nullptr, "Null pointer to mu angle data found in quadrature_data");
  Insist(quad.weights != nullptr, "Null pointer to weight data found in quadrature_data");

  // For 2 and 3-D quadratures, the ordinates have all angles
  if (quad.dimension > 1) {
    Insist(quad.eta != nullptr, "Null pointer to eta angle data found in quadrature_data");
    Insist(quad.xi != nullptr, "Null pointer to xi angle data found in quadrature_data");
  }
}

//------------------------------------------------------------------------------------------------//
void get_quadrature(quadrature_data &quad) {
  using namespace ::rtt_quadrature;

  check_quadrature_validity(quad);

  Ordinate_Set_Factory osf(quad);
  std::shared_ptr<Ordinate_Set> ordinate_set = osf.get_Ordinate_Set();

  vector<Ordinate> const ordinates(ordinate_set->ordinates());
  size_t i(0);
  if (quad.dimension == 1) {
    for (auto ord = ordinates.begin(); ord != ordinates.end(); ++ord, ++i) {
      quad.weights[i] = ord->wt();
      quad.mu[i] = ord->mu();
    }
  } else // dimension == 2 or 3
  {
    for (auto ord = ordinates.begin(); ord != ordinates.end(); ++ord, ++i) {
      quad.weights[i] = ord->wt();
      quad.mu[i] = ord->mu();
      quad.eta[i] = ord->eta();
      quad.xi[i] = ord->xi();
    }
  }

  return;
}

//------------------------------------------------------------------------------------------------//
// end of quadrature/Quadrature_Interface.cc
//------------------------------------------------------------------------------------------------//
