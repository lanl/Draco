//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   quadrature/Quadrature_Interface.cc
 * \author Jae Chang
 * \brief  Quadrature interface definitions
 * \note   Copyright (C) 2015-2026 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "quadrature/Quadrature_Interface.hh"
#include "quadrature/Ordinate.hh"
#include "quadrature/Ordinate_Set.hh"
#include "quadrature/config.h"

#ifdef GSL_FOUND
#include "quadrature/Ordinate_Set_Factory.hh"
#endif

#include "dsxx/Assert.hh"
#include <cstddef> // size_t
#include <memory>
#include <vector>

//! An extern "C" interface to default constructor
void init_quadrature(quadrature_data &quad) { quad = quadrature_data(); }

#ifdef GSL_FOUND
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

  Ordinate_Set_Factory const osf(quad);
  std::shared_ptr<Ordinate_Set> const ordinate_set = osf.get_Ordinate_Set();

  std::vector<Ordinate> const ordinates(ordinate_set->ordinates());
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
#else
// Function to check basic validity of quadrature_data
void check_quadrature_validity(const quadrature_data & /*quad*/) {
  Insist(0, "check_quadrature_validity not supported if GSL is not set");
}

//------------------------------------------------------------------------------------------------//
void get_quadrature(quadrature_data & /*quad*/) {
  Insist(0, "get_quadrature not supported if GSL is not set");
}

#endif
//------------------------------------------------------------------------------------------------//
// end of quadrature/Quadrature_Interface.cc
//------------------------------------------------------------------------------------------------//
