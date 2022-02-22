//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   quadrature/Gauss_Legendre.cc
 * \author Kelly Thompson
 * \date   Tue Feb 22 10:21:50 2000
 * \brief  A class representing an interval Gauss-Legendre quadrature set.
 * \note   Copyright (C) 2012-2022 Triad National Security, LLC., All rights reserved.
 *         reserved.  */
//------------------------------------------------------------------------------------------------//

#include "Gauss_Legendre.hh"
#include "gauleg.hh"
#include "ds++/DracoStrings.hh"
#include <numeric>

namespace rtt_quadrature {
using namespace std;
using rtt_dsxx::to_string;

Gauss_Legendre::Gauss_Legendre(unsigned sn_order) : Interval_Quadrature(sn_order) {
  Require(sn_order > 0 && sn_order % 2 == 0);

  // base class data members
  mu_.resize(sn_order);
  wt_.resize(sn_order);

  double const mu1 = -1; // range of direction
  double const mu2 = 1;
  gauleg(mu1, mu2, mu_, wt_, sn_order_);

  // Sanity Checks: none at present

  Ensure(check_class_invariants());
  Ensure(this->sn_order() == sn_order);
}

//------------------------------------------------------------------------------------------------//
/* virtual */
string Gauss_Legendre::name() const { return "Gauss-Legendre"; }

//------------------------------------------------------------------------------------------------//
/* virtual */
string Gauss_Legendre::parse_name() const { return "gauss legendre"; }

//------------------------------------------------------------------------------------------------//
/* virtual */
unsigned Gauss_Legendre::number_of_levels() const { return sn_order_; }

//------------------------------------------------------------------------------------------------//
/* virtual */ string Gauss_Legendre::as_text(string const &indent) const {
  string Result = indent + "type = gauss legendre" + indent + "  order = " + to_string(sn_order_) +
                  indent + "end";

  return Result;
}

//------------------------------------------------------------------------------------------------//
bool Gauss_Legendre::check_class_invariants() const { return sn_order_ > 0 && sn_order_ % 2 == 0; }

//------------------------------------------------------------------------------------------------//
/* virtual */
vector<Ordinate> Gauss_Legendre::create_level_ordinates_(double const norm) const {
  // Preconditions checked in create_ordinate_set

  unsigned const numPoints(sn_order());

  double sumwt = 0.0;
  for (size_t i = 0; i < numPoints; ++i)
    sumwt += wt_[i];

  double c = norm / sumwt;

  // build the set of ordinates
  vector<Ordinate> Result(numPoints);
  for (size_t i = 0; i < numPoints; ++i) {
    // This is a 1D set.
    Result[i] = Ordinate(mu_[i], c * wt_[i]);
  }

  return Result;
}

} // end namespace rtt_quadrature

//------------------------------------------------------------------------------------------------//
// end of quadrature/Quadrature.cc
//------------------------------------------------------------------------------------------------//
