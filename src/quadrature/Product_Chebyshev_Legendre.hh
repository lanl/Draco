//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   quadrature/Product_Chebyshev_Legendre.hh
 * \author James S. Warsa
 * \brief  A class for Product Chebyshev-Gauss-Legendre quadrature sets.
 * \note   Copyright (C) 2013-2026 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef quadrature_Product_Chebyshev_Legendre_hh
#define quadrature_Product_Chebyshev_Legendre_hh

#include "Octant_Quadrature.hh"

namespace rtt_quadrature {

//================================================================================================//
/*!
 * \class Product_Chebyshev_Legendre
 * \brief A class to encapsulate a triangular Chebyshev-Legendre quadrature set.
 */
//================================================================================================//

class Product_Chebyshev_Legendre : public Octant_Quadrature {
public:
  // CREATORS

  // The default values for snOrder_ and norm_ were set in QuadCreator.
  Product_Chebyshev_Legendre(unsigned sn_order, unsigned azimuthal_order)
      : Octant_Quadrature(sn_order), azimuthal_order_(azimuthal_order) {
    Require(sn_order > 0 && sn_order % 2 == 0);
    Require(azimuthal_order > 0 && azimuthal_order % 2 == 0);
    Ensure(this->azimuthal_order() == azimuthal_order);
  }

  Product_Chebyshev_Legendre(unsigned sn_order, unsigned azimuthal_order, unsigned const mu_axis,
                             unsigned const eta_axis)
      : Octant_Quadrature(sn_order, mu_axis, eta_axis), azimuthal_order_(azimuthal_order) {
    Require(sn_order > 0 && sn_order % 2 == 0);
    Require(azimuthal_order > 0 && azimuthal_order % 2 == 0);

    Ensure(this->azimuthal_order() == azimuthal_order);
  }

  Product_Chebyshev_Legendre(); // disable default construction

  // ACCESSORS

  [[nodiscard]] unsigned azimuthal_order() const {
    Ensure(azimuthal_order_ > 0 && azimuthal_order_ % 2 == 0);
    return azimuthal_order_;
  }

  // SERVICES

  // These functions override the virtual member functions specifed in the parent class Quadrature.

  [[nodiscard]] std::string name() const override;
  [[nodiscard]] std::string parse_name() const override;
  [[nodiscard]] Quadrature_Class quadrature_class() const override;
  [[nodiscard]] unsigned number_of_levels() const override;
  // NOLINTNEXTLINE(misc-override-with-different-visibility)
  [[nodiscard]] std::string as_text(std::string const &indent) const override;

  // STATICS

  static std::shared_ptr<Quadrature> parse(rtt_parser::Token_Stream &tokens);

private:
  // IMPLEMENTATION

  //! Virtual hook for create_ordinate_set
  void create_octant_ordinates_( // NOLINT(misc-override-with-different-visibility)
      std::vector<double> &mu, std::vector<double> &eta, std::vector<double> &wt) const override;

  unsigned const azimuthal_order_;

  // DATA
};

} // end namespace rtt_quadrature

#endif // quadrature_Product_Chebyshev_Legendre_hh

//------------------------------------------------------------------------------------------------//
// end of quadrature/Product_Chebyshev_Legendre.hh
//------------------------------------------------------------------------------------------------//
