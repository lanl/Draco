//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   quadrature/Ordinate_Set.cc
 * \author Kent Budge
 * \date   Tue Dec 21 14:20:03 2004
 * \brief  Declaration file for the class rtt_quadrature::Ordinate.
 * \note   Copyright (C) 2012-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "Ordinate_Set.hh"
#include <algorithm>
#include <iomanip>
#include <iostream>

namespace {
using namespace std;
using namespace rtt_quadrature;
using rtt_dsxx::soft_equiv;

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

  if (soft_equiv(a.eta(), b.eta())) {
    if (soft_equiv(a.mu(), b.mu())) {
      if (soft_equiv(a.xi(), b.xi())) {
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
  else if (!soft_equiv(fabs(a.xi()), fabs(b.xi()), 1.0e-14)) {
    return (fabs(a.xi()) > fabs(b.xi()));
  } else if (!soft_equiv(fabs(a.eta()), fabs(b.eta()), 1.0e-14)) {
    return (fabs(a.eta()) < fabs(b.eta()));
  } else {
    return (!soft_equiv(fabs(a.mu()), fabs(b.mu()), 1.0e-14) && fabs(a.mu()) > fabs(b.mu()));
  }
}

//------------------------------------------------------------------------------------------------//
/*!
 * Construct an Ordinate_Set.
 *
 * \param dimension Dimension of the problem. Must be consistent with the geometry.
 * \param geometry Geometry of the problem.
 * \param ordinates Ordinate set for this problem.
 * \param has_starting_directions Has starting directions on each level set.
 * \param has_extra_starting_directions Has extra directions on each level set. In most geometries,
 *          an additional ordinate is added that is opposite in direction to the starting
 *          direction. This is used to implement reflection exactly in curvilinear coordinates. In
 *          1D spherical, that means an additional angle is added at mu=1. In axisymmetric, that
 *          means additional angles are added that are oriented opposite to the incoming starting
 *          direction on each level.
 * \param ordering Ordering into which to sort the ordinates.
*/
Ordinate_Set::Ordinate_Set(unsigned const dimension, rtt_mesh_element::Geometry const geometry,
                           std::vector<Ordinate> const &ordinates,
                           bool const has_starting_directions,
                           bool const has_extra_starting_directions, Ordering const ordering)
    : geometry_(geometry), dimension_(dimension), has_starting_directions_(has_starting_directions),
      has_extra_starting_directions_(has_extra_starting_directions), ordering_(ordering),
      norm_(0.0), ordinates_(ordinates) {
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

  switch (ordering) {
  case LEVEL_ORDERED:
    sort(ordinates_.begin(), ordinates_.end(), level_compare);
    break;

  case OCTANT_ORDERED:
    sort(ordinates_.begin(), ordinates_.end(), octant_compare);
    break;
  }

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
  using std::endl;
  using std::setprecision;

  if (dimension() == 1 && geometry() != rtt_mesh_element::Geometry::AXISYMMETRIC) {
    cout << endl << "The Quadrature directions and weights are:" << endl << endl;
    cout << "   m  \t    mu        \t     wt      " << endl;
    cout << "  --- \t------------- \t-------------" << endl;
    double sum_wt = 0.0;
    for (size_t ix = 0; ix < ordinates_.size(); ++ix) {
      cout << "   " << setprecision(5) << ix << "\t" << setprecision(10) << ordinates_[ix].mu()
           << "\t" << setprecision(10) << ordinates_[ix].wt() << endl;
      sum_wt += ordinates_[ix].wt();
    }
    cout << endl << "  The sum of the weights is " << sum_wt << endl;
    cout << endl;
  } else {
    cout << endl << "The Quadrature directions and weights are:" << endl << endl;
    cout << "   m  \t    mu        \t    eta       \t    xi        \t     wt   "
            "   "
         << endl;
    cout << "  --- \t------------- \t------------- \t------------- "
            "\t-------------"
         << endl;
    double sum_wt = 0.0;
    for (size_t ix = 0; ix < ordinates_.size(); ++ix) {
      cout << "   " << ix << "\t" << setprecision(10) << ordinates_[ix].mu() << "\t"
           << setprecision(10) << ordinates_[ix].eta() << "\t" << setprecision(10)
           << ordinates_[ix].xi() << "\t" << setprecision(10) << ordinates_[ix].wt() << endl;
      sum_wt += ordinates_[ix].wt();
    }
    cout << endl << "  The sum of the weights is " << sum_wt << endl;
    cout << endl;
  }
}

} // end namespace rtt_quadrature

//------------------------------------------------------------------------------------------------//
// end of quadrature/Ordinate_Set.cc
//------------------------------------------------------------------------------------------------//
