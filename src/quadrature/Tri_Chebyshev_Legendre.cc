//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   quadrature/Tri_Chebyshev_Legendre.cc
 * \author Kelly Thompson
 * \note   Copyright (C) 2012-2026 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "quadrature/Tri_Chebyshev_Legendre.hh"
#include "dsxx/DracoStrings.hh"
#include "quadrature/Gauss_Legendre.hh"
#include "quadrature/Octant_Quadrature.hh"
#include "quadrature/Quadrature_Class.hh"
#include <cmath>   // sin, cos, sqrt
#include <cstddef> // size_t
#include <memory>  // shared_ptr
#include <numbers>
#include <string>
#include <vector>

namespace rtt_quadrature {
using namespace rtt_dsxx;

//------------------------------------------------------------------------------------------------//
std::string Tri_Chebyshev_Legendre::name() const { return "Tri Chebyshev Legendre"; }

//------------------------------------------------------------------------------------------------//
std::string Tri_Chebyshev_Legendre::parse_name() const { return "tri cl"; }

//------------------------------------------------------------------------------------------------//
Quadrature_Class Tri_Chebyshev_Legendre::quadrature_class() const { return TRIANGLE_QUADRATURE; }

//------------------------------------------------------------------------------------------------//
unsigned Tri_Chebyshev_Legendre::number_of_levels() const { return sn_order_; }

//------------------------------------------------------------------------------------------------//
std::string Tri_Chebyshev_Legendre::as_text(std::string const &indent) const {
  return indent + "type = tri cl" + indent + "  order = " + to_string(sn_order_) +
         Octant_Quadrature::as_text(indent);
}

//------------------------------------------------------------------------------------------------//
void Tri_Chebyshev_Legendre::create_octant_ordinates_(std::vector<double> &mu,
                                                      std::vector<double> &eta,
                                                      std::vector<double> &wt) const {

  // The number of quadrature levels is equal to the requested SN order.
  size_t const levels = sn_order_;

  // We build the 3-D first, then edit as appropriate.
  size_t const numOrdinates = levels * (levels + 2) / 8;

  // Force the direction vectors to be the correct length.
  mu.resize(numOrdinates);
  eta.resize(numOrdinates);
  wt.resize(numOrdinates);

  std::shared_ptr<Gauss_Legendre> const GL(new Gauss_Legendre(sn_order_));

  // NOTE: this aligns the gauss points with the x-axis (r-axis in cylindrical coords)
  unsigned icount = 0;

  for (unsigned i = 0; i < levels / 2; ++i) {
    double const xmu = GL->mu(i);
    double const xwt = GL->wt(i);
    double const xsr = sqrt(1.0 - xmu * xmu);

    unsigned const k = 2 * (i + 1);

    for (unsigned j = 0; j < k / 2; ++j) {
      unsigned const ordinate = icount;

      mu[ordinate] = xsr * cos(std::numbers::pi * (2.0 * j + 1.0) / k / 2.0);
      eta[ordinate] = xsr * sin(std::numbers::pi * (2.0 * j + 1.0) / k / 2.0);
      wt[ordinate] = xwt / k;

      ++icount;
    }
  }
}

} // end namespace rtt_quadrature

//------------------------------------------------------------------------------------------------//
// end of Tri_Chebyshev_Legendre.cc
//------------------------------------------------------------------------------------------------//
