//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   quadrature/Product_Chebyshev_Legendre.cc
 * \author James S. Warsa
 * \brief  A class for Product Chebyshev-Gauss-Legendre quadrature sets.
 * \note   Copyright (C) 2013-2026 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "quadrature/Product_Chebyshev_Legendre.hh"
#include "dsxx/DracoStrings.hh"
#include "quadrature/Gauss_Legendre.hh"
#include "quadrature/Octant_Quadrature.hh"
#include "quadrature/Quadrature_Class.hh"
#include <cmath>   // sin
#include <cstddef> // size_t
#include <memory>  // shared_ptr
#include <numbers>
#include <string>
#include <vector>

namespace rtt_quadrature {
using namespace rtt_dsxx;

//------------------------------------------------------------------------------------------------//
std::string Product_Chebyshev_Legendre::name() const { return "Product Chebyshev Legendre"; }

//------------------------------------------------------------------------------------------------//
std::string Product_Chebyshev_Legendre::parse_name() const { return "product cl"; }

//------------------------------------------------------------------------------------------------//
Quadrature_Class Product_Chebyshev_Legendre::quadrature_class() const { return OCTANT_QUADRATURE; }

//------------------------------------------------------------------------------------------------//
unsigned Product_Chebyshev_Legendre::number_of_levels() const { return sn_order_; }

//------------------------------------------------------------------------------------------------//
std::string Product_Chebyshev_Legendre::as_text(std::string const &indent) const {
  return indent + "type = " + parse_name() + indent + "  order = " + to_string(sn_order()) + " " +
         to_string(azimuthal_order_) + Octant_Quadrature::as_text(indent);
}

//------------------------------------------------------------------------------------------------//
void Product_Chebyshev_Legendre::create_octant_ordinates_(std::vector<double> &mu,
                                                          std::vector<double> &eta,
                                                          std::vector<double> &wt) const {
  using rtt_dsxx::soft_equiv;
  using std::cos;
  using std::fabs;
  using std::sqrt;

  // The number of quadrature levels is equal to the requested SN order.
  size_t const levels = sn_order();

  // We build the 3-D first, then edit as appropriate.

  size_t const numOrdinates = levels * azimuthal_order_ / 4;

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

    for (unsigned j = 0; j < azimuthal_order_ / 2; ++j) {
      unsigned const ordinate = icount;

      mu[ordinate] = xsr * cos(std::numbers::pi * (2.0 * j + 1.0) / azimuthal_order_ / 2.0);
      eta[ordinate] = xsr * sin(std::numbers::pi * (2.0 * j + 1.0) / azimuthal_order_ / 2.0);
      wt[ordinate] = xwt / azimuthal_order_;

      ++icount;
    }
  }
}

} // end namespace rtt_quadrature

//------------------------------------------------------------------------------------------------//
// end of Product_Chebyshev_Legendre.cc
//------------------------------------------------------------------------------------------------//
