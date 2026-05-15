//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   quadrature/Ordinate_Set.cc
 * \author Kent Budge
 * \brief  Declaration file for the class rtt_quadrature::Ordinate.
 * \note   Copyright (C) 2012-2026 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "quadrature/Ordinate_Set.hh"
#include "dsxx/Assert.hh"
#include "dsxx/Soft_Equivalence.hh"
#include "mesh_element/Geometry.hh"
#include "quadrature/Ordinate.hh"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <vector>

namespace {
using namespace std;
using namespace rtt_quadrature;

// convenience functions to check ordinates

//------------------------------------------------------------------------------------------------//
bool check_4(std::vector<Ordinate> const &ordinates) {
  // In 1-D spherical geometry, the ordinates must be confined to the first two octants.

  size_t const N = ordinates.size();
  for (unsigned i = 0; i < N; ++i) {
    if (ordinates[i].eta() < 0 || ordinates[i].xi() < 0)
      return false;
  }
  return true;
}

//------------------------------------------------------------------------------------------------//
bool check_2(std::vector<Ordinate> const &ordinates) {
  // In 2-D geometry, the ordinates must be confined to the first four octants

  size_t const N = ordinates.size();
  for (unsigned i = 0; i < N; ++i) {
    if (ordinates[i].xi() < 0)
      return false;
  }
  return true;
}

} // end anonymous namespace

namespace rtt_quadrature {

//------------------------------------------------------------------------------------------------//
/* static */
bool Ordinate_Set::level_compare(Ordinate const &a, Ordinate const &b) {
  // Note that x==r==mu, z==xi

  if (rtt_dsxx::soft_equiv(a.eta(), b.eta())) {
    if (rtt_dsxx::soft_equiv(a.mu(), b.mu())) {
      if (rtt_dsxx::soft_equiv(a.xi(), b.xi())) {
        return false;
      } else {
        return a.xi() < b.xi();
      }
    } else {
      return a.mu() < b.mu();
    }
  } else {
    return a.eta() < b.eta();
  }
}

//------------------------------------------------------------------------------------------------//
bool Ordinate_Set::octant_compare(Ordinate const &a, Ordinate const &b) {
  // We initially sort by octant. Only the +++ octant is actually used by PARTISN-type sweepers that
  // assume all quadratures are octant quadratures.

  if (a.xi() < 0 && b.xi() > 0) {
    return true;
  } else if (a.xi() > 0 && b.xi() < 0) {
    return false;
  } else if (a.eta() < 0 && b.eta() > 0) {
    return true;
  } else if (a.eta() > 0 && b.eta() < 0) {
    return false;
  } else if (a.mu() < 0 && b.mu() > 0) {
    return true;
  } else if (a.mu() > 0 && b.mu() < 0) {
    return false;
  }
  // Within an octant, we sort by decreasing absolute xi, then increasing absolute eta, to be
  // consistent with PARTISN.
  else if (!rtt_dsxx::soft_equiv(fabs(a.xi()), fabs(b.xi()), 1.0e-14)) {
    return (fabs(a.xi()) > fabs(b.xi()));
  } else if (!rtt_dsxx::soft_equiv(fabs(a.eta()), fabs(b.eta()), 1.0e-14)) {
    return (fabs(a.eta()) < fabs(b.eta()));
  } else {
    return (!rtt_dsxx::soft_equiv(fabs(a.mu()), fabs(b.mu()), 1.0e-14) &&
            fabs(a.mu()) > fabs(b.mu()));
  }
}

//------------------------------------------------------------------------------------------------//
/*!
 * Construct an Ordinate_Set.
 *
 * \param[in] dimension Dimension of the problem. Must be consistent with the geometry.
 * \param[in] geometry Geometry of the problem.
 * \param ordinates Ordinate set for this problem.
 * \param[in] has_starting_directions Has starting directions on each level set.
 * \param[in] has_extra_starting_directions Has extra directions on each level set. In most
 *          geometries, an additional ordinate is added that is opposite in direction to the
 *          starting direction. This is used to implement reflection exactly in curvilinear
 *          coordinates. In 1D spherical, that means an additional angle is added at mu=1. In
 *          axisymmetric, that means additional angles are added that are oriented opposite to the
 *          incoming starting direction on each level.
 * \param[in] ordering Ordering into which to sort the ordinates.
*/
Ordinate_Set::Ordinate_Set(unsigned const dimension, rtt_mesh_element::Geometry const geometry,
                           std::vector<Ordinate> const &ordinates,
                           bool const has_starting_directions,
                           bool const has_extra_starting_directions, Ordering const ordering)
    : geometry_(geometry), dimension_(dimension), has_starting_directions_(has_starting_directions),
      has_extra_starting_directions_(has_extra_starting_directions), ordering_(ordering),
      ordinates_(ordinates) {
  Require(dimension >= 1 && dimension <= 3);
  Require(geometry != rtt_mesh_element::Geometry::AXISYMMETRIC || dimension < 3);
  Require(geometry != rtt_mesh_element::Geometry::SPHERICAL || dimension < 2);
  Require(has_starting_directions || !has_extra_starting_directions);
  Require(dimension > 1 || geometry == rtt_mesh_element::Geometry::SPHERICAL || check_4(ordinates));
  Require(dimension != 2 || check_2(ordinates));

  norm_ = 0.0;
  size_t const N = ordinates_.size();
  for (unsigned i = 0; i < N; ++i) {
    norm_ += ordinates[i].wt();
  }

#if defined(__clang__)
#if defined(__clang_major__) && __clang_major__ > 17
// llvm warns about the missing 'default:' label, but it also warns if this label is present but
// can never be reached.
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wswitch-default"
#endif
#endif
  switch (ordering) {
  case LEVEL_ORDERED:
    std::ranges::sort(ordinates_, level_compare);
    break;

  case OCTANT_ORDERED:
    std::ranges::sort(ordinates_, octant_compare);
    break;
  }
#if defined(__clang__)
#if defined(__clang_major__) && __clang_major__ > 17
#pragma clang diagnostic pop
#endif
#endif

  Ensure(check_class_invariants());
  Ensure(this->has_starting_directions() == has_starting_directions);
  Ensure(this->has_extra_starting_directions() == has_extra_starting_directions);
  Ensure(this->ordering() == ordering);
}

//! Copy constructor
Ordinate_Set::Ordinate_Set(Ordinate_Set const &other)
    : geometry_(other.geometry()), dimension_(other.dimension()),
      has_starting_directions_(other.has_starting_directions()),
      has_extra_starting_directions_(other.has_extra_starting_directions()),
      ordering_(other.ordering()), norm_(other.norm()), ordinates_(other.ordinates()) {}

//------------------------------------------------------------------------------------------------//
bool Ordinate_Set::check_class_invariants() const {
  return (dimension_ >= 1 && dimension_ <= 3) &&
         (geometry_ != rtt_mesh_element::Geometry::AXISYMMETRIC || dimension_ < 3) &&
         (geometry_ != rtt_mesh_element::Geometry::SPHERICAL || dimension_ < 2) &&
         (has_starting_directions_ || !has_extra_starting_directions_) &&
         (dimension_ > 1 || geometry_ == rtt_mesh_element::Geometry::SPHERICAL ||
          check_4(ordinates_)) &&
         (dimension_ != 2 || check_2(ordinates_));
}

//------------------------------------------------------------------------------------------------//
void Ordinate_Set::display() const {
  using std::cout;
  using std::setprecision;

  if (dimension() == 1 && geometry() != rtt_mesh_element::Geometry::AXISYMMETRIC) {
    cout << "\nThe Quadrature directions and weights are:\n\n"
         << "   m  \t    mu        \t     wt      \n"
         << "  --- \t------------- \t-------------\n";
    double sum_wt = 0.0;
    for (size_t ix = 0; ix < ordinates_.size(); ++ix) {
      cout << "   " << setprecision(5) << ix << "\t" << setprecision(10) << ordinates_[ix].mu()
           << "\t" << setprecision(10) << ordinates_[ix].wt() << '\n';
      sum_wt += ordinates_[ix].wt();
    }
    cout << "\n  The sum of the weights is " << sum_wt << "\n\n";
  } else {
    cout << "\nThe Quadrature directions and weights are:\n\n"
         << "   m  \t    mu        \t    eta       \t    xi        \t     wt      \n"
         << "  --- \t------------- \t------------- \t------------- \t-------------\n";
    double sum_wt = 0.0;
    for (size_t ix = 0; ix < ordinates_.size(); ++ix) {
      cout << "   " << ix << "\t" << setprecision(10) << ordinates_[ix].mu() << "\t"
           << setprecision(10) << ordinates_[ix].eta() << "\t" << setprecision(10)
           << ordinates_[ix].xi() << "\t" << setprecision(10) << ordinates_[ix].wt() << '\n';
      sum_wt += ordinates_[ix].wt();
    }
    cout << "\n  The sum of the weights is " << sum_wt << "\n\n";
  }
}

} // end namespace rtt_quadrature

//------------------------------------------------------------------------------------------------//
// end of quadrature/Ordinate_Set.cc
//------------------------------------------------------------------------------------------------//
