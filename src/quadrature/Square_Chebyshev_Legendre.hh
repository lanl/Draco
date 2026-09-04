//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   quadrature/Square_Chebyshev_Legendre.hh
 * \author Kelly Thompson
 * \brief  A class to encapsulate a 3D Level Symmetric quadrature set.
 * \note   Copyright (C) 2012-2026 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef quadrature_Square_Chebyshev_Legendre_hh
#define quadrature_Square_Chebyshev_Legendre_hh

#include "Octant_Quadrature.hh"

namespace rtt_quadrature {

//================================================================================================//
/*!
 * \class Square_Chebyshev_Legendre
 * \brief A class to encapsulate a triangular Chebyshev-Legendre quadrature set.
 */
//================================================================================================//

class Square_Chebyshev_Legendre : public Octant_Quadrature {
public:
  // CREATORS

  // The default values for snOrder_ and norm_ were set in QuadCreator.
  explicit Square_Chebyshev_Legendre(unsigned sn_order)
      : Octant_Quadrature(sn_order)

  {
    Require(sn_order > 0 && sn_order % 2 == 0);
  }

  Square_Chebyshev_Legendre(unsigned sn_order, unsigned const mu_axis, unsigned const eta_axis)
      : Octant_Quadrature(sn_order, mu_axis, eta_axis)

  {
    Require(sn_order > 0 && sn_order % 2 == 0);
  }

  Square_Chebyshev_Legendre(); // disable default construction

  // ACCESSORS

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
};

} // end namespace rtt_quadrature

#endif // quadrature_Square_Chebyshev_Legendre_hh

//------------------------------------------------------------------------------------------------//
// end of quadrature/Square_Chebyshev_Legendre.hh
//------------------------------------------------------------------------------------------------//
