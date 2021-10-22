//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   quadrature/Octant_Quadrature.cc
 * \author Kent Budge
 * \date   Friday, Nov 30, 2012, 08:27 am
 * \brief  Implementation for Octant_Quadrature
 * \note   Copyright (C) 2012-2021 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "Octant_Quadrature.hh"
#include "ds++/DracoStrings.hh"
#include "ds++/Soft_Equivalence.hh"
#include "units/PhysicalConstants.hh"
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>

namespace rtt_quadrature {
using namespace rtt_dsxx;

//------------------------------------------------------------------------------------------------//

bool Octant_Quadrature::has_axis_assignments() const { return has_axis_assignments_; }

//------------------------------------------------------------------------------------------------//
std::vector<Ordinate> Octant_Quadrature::create_ordinates_(
    unsigned const dimension, Geometry const geometry, double const norm, unsigned const mu_axis,
    unsigned const eta_axis, bool const include_starting_directions,
    bool const include_extra_directions) const {
  using rtt_dsxx::soft_equiv;

  // We build the 3-D first, then edit as appropriate.

  std::vector<double> mu, eta, wt;

  create_octant_ordinates_(mu, eta, wt);

  size_t const octantOrdinates = mu.size();
  Check(octantOrdinates > 0);
  Check(eta.size() == octantOrdinates);
  Check(wt.size() == octantOrdinates);

  size_t numOrdinates = octantOrdinates * 8;
  mu.resize(numOrdinates);
  eta.resize(numOrdinates);
  wt.resize(numOrdinates);

  // Evaluate mu and eta for octants 2-4
  for (size_t octant = 2; octant <= 4; ++octant)
    for (size_t n = 0; n <= octantOrdinates - 1; ++n) {
      size_t const m = (octant - 1) * octantOrdinates + n;
      Check(m < mu.size() && m < eta.size() && m < wt.size());
      switch (octant) {
      case 2:
        mu[m] = -mu[n];
        eta[m] = eta[n];
        wt[m] = wt[n];
        break;

      case 3:
        mu[m] = -mu[n];
        eta[m] = -eta[n];
        wt[m] = wt[n];
        break;

      case 4:
        mu[m] = mu[n];
        eta[m] = -eta[n];
        wt[m] = wt[n];
        break;
      default:
        Insist(false, "Octant value should only be 2, 3 or 4 in this loop.");
        break;
      }
    }

  // Evaluate mu and eta for octants 5-8
  for (size_t n = 0; n <= 4 * octantOrdinates - 1; ++n) {
    Check(n + 4 * octantOrdinates < mu.size() && n + 4 * octantOrdinates < eta.size() &&
          n + 4 * octantOrdinates < wt.size());
    mu[n + 4 * octantOrdinates] = mu[n];
    eta[n + 4 * octantOrdinates] = eta[n];
    wt[n + 4 * octantOrdinates] = wt[n];
  }

  // Evaluate xi for all octants
  std::vector<double> xi(numOrdinates);
  for (size_t n = 0; n <= 4 * octantOrdinates - 1; ++n) {
    Check(n < xi.size());
    xi[n] = std::sqrt(1.0 - (mu[n] * mu[n] + eta[n] * eta[n]));
  }

  for (size_t n = 0; n <= 4 * octantOrdinates - 1; ++n) {
    Check(n + 4 * octantOrdinates < xi.size());
    xi[n + 4 * octantOrdinates] = -xi[n];
  }

  std::vector<Ordinate> Result;

  if (dimension == 3) {
    map_axes_(mu_axis, eta_axis, mu, eta, xi);

    // Copy all ordinates into the result vector.
    Result.resize(numOrdinates);
    for (size_t i = 0; i < numOrdinates; ++i) {
      Result[i] = Ordinate(mu[i], eta[i], xi[i], wt[i]);
    }
  } else if (dimension == 2 || geometry == rtt_mesh_element::Geometry::AXISYMMETRIC) {
    map_axes_(mu_axis, eta_axis, mu, eta, xi);

    // Copy the half-sphere
    Result.resize(4 * octantOrdinates);
    unsigned m = 0;
    for (size_t i = 0; i < numOrdinates; ++i) {
      if (xi[i] > 0.0 && (dimension > 1 || eta[i] > 0.0)) {
        Check(m < Result.size());
        Result[m++] = Ordinate(mu[i], eta[i], xi[i], wt[i]);
      }
    }
    Result.resize(m);

    // Add starting directions if appropriate
    add_2D_starting_directions_(geometry, include_starting_directions, include_extra_directions,
                                Result);
  } else {
    Check(dimension == 1 && geometry != rtt_mesh_element::Geometry::AXISYMMETRIC);
    Check(mu_axis == 2);

    map_axes_(0, 2, mu, eta, xi);

    // Only need the quarter sphere
    Result.resize(numOrdinates / 4);
    unsigned m = 0;
    for (size_t i = 0; i < numOrdinates; ++i) {
      if (mu[i] > 0.0 && xi[i] > 0.0) {
        Check(m < Result.size());
        Result[m++] = Ordinate(mu[i], eta[i], xi[i], wt[i]);
      }
    }

    // Sort

    std::sort(Result.begin(), Result.end(), Ordinate_Set::level_compare);

    // Now sum around the axis.
    m = 0;
    numOrdinates /= 4;
    double eta0 = Result[0].eta();
    double sum = Result[0].wt();
    for (unsigned i = 1; i < numOrdinates; ++i) {
      double old_eta = eta0;
      eta0 = Result[i].eta();
      if (!soft_equiv(eta0, old_eta)) {
        // New level
        Result[m++] = Ordinate(old_eta, sum);
        sum = Result[i].wt();
      } else {
        // Still on old level
        sum += Result[i].wt();
      }
    }
    // Final level
    Result[m++] = Ordinate(eta0, sum);
    numOrdinates = m;
    Result.resize(numOrdinates);

    // Add starting directions if appropriate
    add_1D_starting_directions_(geometry, include_starting_directions, include_extra_directions,
                                Result);
  }

  numOrdinates = Result.size();

  // Normalize the quadrature set
  double wsum = 0.0;
  for (size_t n = 0; n <= numOrdinates - 1; ++n)
    wsum = wsum + Result[n].wt();

  if (dimension == 1 && geometry != rtt_mesh_element::Geometry::AXISYMMETRIC) {
    for (size_t n = 0; n <= numOrdinates - 1; ++n)
      Result[n] = Ordinate(Result[n].mu(), Result[n].wt() * (norm / wsum));
  } else {
    for (size_t n = 0; n <= numOrdinates - 1; ++n)
      Result[n] =
          Ordinate(Result[n].mu(), Result[n].eta(), Result[n].xi(), Result[n].wt() * (norm / wsum));
  }
  return Result;
}

//------------------------------------------------------------------------------------------------//
std::vector<Ordinate> Octant_Quadrature::create_ordinates_(unsigned dimension, Geometry geometry,
                                                           double norm,
                                                           bool include_starting_directions,
                                                           bool include_extra_directions) const {
  unsigned mu_axis(0), eta_axis(0);
  if (has_axis_assignments_) {
    mu_axis = mu_axis_;
    eta_axis = eta_axis_;
  } else {
    switch (dimension) {
    case 1:
      switch (geometry) {
      case rtt_mesh_element::Geometry::AXISYMMETRIC:
        mu_axis = 0;
        eta_axis = 2;
        break;

      default:
        mu_axis = 2;
        eta_axis = 1;
        break;
      }
      break;

    case 2:
      switch (geometry) {
      case rtt_mesh_element::Geometry::AXISYMMETRIC:
        mu_axis = 0;
        eta_axis = 2;
        break;

      default:
        mu_axis = 0;
        eta_axis = 1;
        break;
      }
      break;

    case 3:
      mu_axis = 0;
      eta_axis = 1;
      break;

    default:
      Insist(false, "bad case");
    }
  }
  return create_ordinates_(dimension, geometry, norm, mu_axis, eta_axis,
                           include_starting_directions, include_extra_directions);
}

//------------------------------------------------------------------------------------------------//
/*!
 * Pure virtual used in conjuction with child implementations, for common features.
 */
std::string Octant_Quadrature::as_text(std::string const &indent) const {
  std::string Result;

  if (has_axis_assignments_) {
    Result += indent + "  axis assignments, mu = " + to_string(mu_axis_) +
              " eta = " + to_string(eta_axis_);
  }

  Result += indent + "end";

  return Result;
}

} // end namespace rtt_quadrature

//------------------------------------------------------------------------------------------------//
// end of Octant_Quadrature.cc
//------------------------------------------------------------------------------------------------//
