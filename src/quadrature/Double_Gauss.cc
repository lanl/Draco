//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   quadrature/Double_Gauss.cc
 * \author Kelly Thompson
 * \brief  A class representing an interval double Gauss-Legendre quadrature set.
 * \note   Copyright (C) 2012-2026 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "quadrature/Double_Gauss.hh"
#include "dsxx/Assert.hh"
#include "dsxx/DracoStrings.hh"
#include "quadrature/Gauss_Legendre.hh"
#include "quadrature/Interval_Quadrature.hh"
#include "quadrature/Ordinate.hh"
#include <cstddef> // size_t
#include <memory>  // shared_ptr
#include <string>
#include <vector>

namespace rtt_quadrature {
using namespace std;
using rtt_dsxx::to_string;

Double_Gauss::Double_Gauss(unsigned sn_order) : Interval_Quadrature(sn_order) {
  Require(sn_order > 0 && sn_order % 2 == 0);

  // base class data members
  mu_.resize(sn_order);
  wt_.resize(sn_order);

  unsigned const numGaussPoints = sn_order;
  unsigned const n(numGaussPoints);
  unsigned const n2(n / 2);

  if (n2 == 1) // 2-point double Gauss is just Gauss
  {
    Check(sn_order == 2);

    std::shared_ptr<Gauss_Legendre> const GL(new Gauss_Legendre(sn_order));
    for (unsigned m = 0; m < sn_order; ++m) {
      mu_[m] = GL->mu(m);
      wt_[m] = GL->wt(m);
    }
  } else // Create an N/2-point Gauss quadrature on [-1,1]
  {

    Check(n2 % 2 == 0);
    Check(n2 > 2);

    std::shared_ptr<Gauss_Legendre> const GL(new Gauss_Legendre(n2));

    // map the quadrature onto the two half-ranges

    for (unsigned m = 0; m < n2; ++m) {
      // Map onto [-1,0] then skew-symmetrize (ensuring ascending order on [-1, 1])

      mu_[m] = 0.5 * (GL->mu(m) - 1.0);
      wt_[m] = 0.5 * GL->wt(m);

      mu_[n - m - 1] = -mu_[m];
      wt_[n - m - 1] = wt_[m];
    }
  }

  Ensure(check_class_invariants());
  Ensure(this->sn_order() == sn_order);
}

//------------------------------------------------------------------------------------------------//
/* virtual */
string Double_Gauss::name() const { return "Double-Gauss"; }

//------------------------------------------------------------------------------------------------//
/* virtual */
string Double_Gauss::parse_name() const { return "double gauss"; }

//------------------------------------------------------------------------------------------------//
/* virtual */
unsigned Double_Gauss::number_of_levels() const { return sn_order(); }

//------------------------------------------------------------------------------------------------//
/* virtual */ string Double_Gauss::as_text(string const &indent) const {
  string Result = indent + "type = double gauss" + indent + "  order = " + to_string(sn_order()) +
                  indent + "end";

  return Result;
}

//------------------------------------------------------------------------------------------------//
bool Double_Gauss::check_class_invariants() const { return sn_order() > 0 && sn_order() % 2 == 0; }

//------------------------------------------------------------------------------------------------//
/* virtual */
vector<Ordinate> Double_Gauss::create_level_ordinates_(double const norm) const {
  // Preconditions checked in create_ordinate_set

  unsigned const numPoints(sn_order());

  double sumwt = 0.0;
  for (size_t i = 0; i < numPoints; ++i)
    sumwt += wt_[i];

  double const c = norm / sumwt;

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
// end of quadrature/Double_Gauss.cc
//------------------------------------------------------------------------------------------------//
