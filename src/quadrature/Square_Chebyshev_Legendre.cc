//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   quadrature/Square_Chebyshev_Legendre.cc
 * \author Kelly Thompson
 * \date   Wed Sep  1 10:19:52 2004
 * \note   Copyright (C) 2012-2022 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

#include "Square_Chebyshev_Legendre.hh"
#include "Gauss_Legendre.hh"
#include "ds++/DracoStrings.hh"
#include "units/MathConstants.hh"

namespace rtt_quadrature {
using namespace rtt_dsxx;

//------------------------------------------------------------------------------------------------//
std::string Square_Chebyshev_Legendre::name() const { return "Square Chebyshev Legendre"; }

//------------------------------------------------------------------------------------------------//
std::string Square_Chebyshev_Legendre::parse_name() const { return "square cl"; }

//------------------------------------------------------------------------------------------------//
Quadrature_Class Square_Chebyshev_Legendre::quadrature_class() const { return SQUARE_QUADRATURE; }

//------------------------------------------------------------------------------------------------//
unsigned Square_Chebyshev_Legendre::number_of_levels() const { return sn_order_; }

//------------------------------------------------------------------------------------------------//
std::string Square_Chebyshev_Legendre::as_text(std::string const &indent) const {
  return indent + "type = square cl" + indent + "  order = " + to_string(sn_order()) +
         Octant_Quadrature::as_text(indent);
}

//------------------------------------------------------------------------------------------------//
void Square_Chebyshev_Legendre::create_octant_ordinates_(std::vector<double> &mu,
                                                         std::vector<double> &eta,
                                                         std::vector<double> &wt) const {
  using rtt_dsxx::soft_equiv;
  using std::cos;
  using std::fabs;
  using std::sqrt;

  // The number of quadrature levels is equal to the requested SN order.
  size_t levels = sn_order();

  // We build the 3-D first, then edit as appropriate.

  size_t numOrdinates = levels * levels / 4;

  // Force the direction vectors to be the correct length.
  mu.resize(numOrdinates);
  eta.resize(numOrdinates);
  wt.resize(numOrdinates);

  std::shared_ptr<Gauss_Legendre> GL(new Gauss_Legendre(sn_order_));

  // NOTE: this aligns the gauss points with the x-axis (r-axis in cylindrical coords)

  for (unsigned i = 0; i < levels / 2; ++i) {
    double xmu = GL->mu(i);
    double xwt = GL->wt(i);
    double xsr = sqrt(1.0 - xmu * xmu);

    for (unsigned j = 0; j < levels / 2; ++j) {
      size_t ordinate = j + i * levels / 2;

      mu[ordinate] = xsr * cos(rtt_units::PI * static_cast<double>(2 * j + 1) /
                               static_cast<double>(levels) / 2.0);
      eta[ordinate] = xsr * sin(rtt_units::PI * static_cast<double>(2 * j + 1) /
                                static_cast<double>(levels) / 2.0);
      wt[ordinate] = xwt / static_cast<double>(levels);
    }
  }
}

} // end namespace rtt_quadrature

//------------------------------------------------------------------------------------------------//
// end of Square_Chebyshev_Legendre.cc
//------------------------------------------------------------------------------------------------//
