//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   quadrature/General_Octant_Quadrature.cc
 * \author Kelly Thompson
 * \date   Wed Sep  1 10:19:52 2004
 * \note   Copyright (C) 2012-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "General_Octant_Quadrature.hh"
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
General_Octant_Quadrature::General_Octant_Quadrature(
    unsigned const sn_order, std::vector<double> const &mu, std::vector<double> const &eta,
    std::vector<double> const &xi, std::vector<double> const &wt, unsigned number_of_levels,
    Quadrature_Class const quadrature_class)
    : Octant_Quadrature(sn_order), mu_(mu), eta_(eta), xi_(xi), wt_(wt),
      number_of_levels_(number_of_levels), quadrature_class_(quadrature_class),
      is_open_interval_(soft_equiv(*std::min_element(xi.begin(), xi.end()), -1.0)) {
  Require(mu.size() > 0 && eta.size() == mu.size() && xi.size() == mu.size() &&
          wt.size() == mu.size());

  Require(quadrature_class != INTERVAL_QUADRATURE);

  Require(quadrature_class != TRIANGLE_QUADRATURE ||
          number_of_levels * (number_of_levels + 2) == 8 * mu.size());

  Require(quadrature_class != SQUARE_QUADRATURE ||
          2 * number_of_levels * number_of_levels == 8 * mu.size());

  // is_open_interval_ =
  //     soft_equiv(*std::min_element(xi.begin(), xi.end()), -1.0);

  Ensure(check_class_invariants());
  Ensure(this->mu() == mu);
  Ensure(this->eta() == eta);
  Ensure(this->xi() == xi);
  Ensure(this->wt() == wt);
  //! \bug [clang-tidy] Calling a pure virtual method during construction or destruction is
  //       undefined behavior.  Ensure(this->number_of_levels() == number_of_levels);
  // Ensure(this->quadrature_class() == quadrature_class);
}

//------------------------------------------------------------------------------------------------//
bool General_Octant_Quadrature::check_class_invariants() const {
  return (mu_.size() > 0 && eta_.size() == mu_.size() && xi_.size() == mu_.size() &&
          wt_.size() == mu_.size()) &&

         (quadrature_class_ != INTERVAL_QUADRATURE) &&

         (quadrature_class_ != TRIANGLE_QUADRATURE ||
          number_of_levels_ * (number_of_levels_ + 2) == 8 * mu_.size()) &&

         (quadrature_class_ != SQUARE_QUADRATURE ||
          2 * number_of_levels_ * number_of_levels_ == 8 * mu_.size());
}

//------------------------------------------------------------------------------------------------//
std::string General_Octant_Quadrature::name() const { return "General Octant Quadrature"; }

//------------------------------------------------------------------------------------------------//
std::string General_Octant_Quadrature::parse_name() const { return "general octant quadrature"; }

//------------------------------------------------------------------------------------------------//
Quadrature_Class General_Octant_Quadrature::quadrature_class() const { return quadrature_class_; }

//------------------------------------------------------------------------------------------------//
unsigned General_Octant_Quadrature::number_of_levels() const { return number_of_levels_; }

//------------------------------------------------------------------------------------------------//
std::string General_Octant_Quadrature::as_text(std::string const &indent) const {
  std::string Result = indent + "  type = general octant quadrature";
  Result += indent + "  sn order = " + to_string(sn_order());
  Result += indent + "  number of ordinates = " + to_string(mu_.size());
  Result += indent + "  number of levels = " + to_string(number_of_levels_);

  size_t const N = mu_.size();
  for (size_t i = 0; i < N; ++i) {
    Result += indent + "  " + to_string(mu_[i]);
    Result += "  " + to_string(eta_[i]);
    Result += "  " + to_string(xi_[i]);
    Result += "  " + to_string(wt_[i]);
  }

  Result += indent + "end";

  return Result;
}

//------------------------------------------------------------------------------------------------//
void General_Octant_Quadrature::create_octant_ordinates_(std::vector<double> &mu,
                                                         std::vector<double> &eta,
                                                         std::vector<double> &wt) const {
  mu = mu_;
  eta = eta_;
  wt = wt_;
}

//------------------------------------------------------------------------------------------------//
bool General_Octant_Quadrature::is_open_interval() const { return is_open_interval_; }

} // end namespace rtt_quadrature

//------------------------------------------------------------------------------------------------//
// end of General_Octant_Quadrature.cc
//------------------------------------------------------------------------------------------------//
