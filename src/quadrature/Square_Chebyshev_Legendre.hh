//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   quadrature/Square_Chebyshev_Legendre.hh
 * \author Kelly Thompson
 * \date   Wed Sep  1 10:19:52 2004
 * \brief  A class to encapsulate a 3D Level Symmetric quadrature set.
 * \note   Copyright (C) 2012-2022 Triad National Security, LLC., All rights reserved. */
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

  std::string name() const override;
  std::string parse_name() const override;
  Quadrature_Class quadrature_class() const override;
  unsigned number_of_levels() const override;
  std::string as_text(std::string const &indent) const override;

  // STATICS

  static std::shared_ptr<Quadrature> parse(Token_Stream &tokens);

private:
  // IMPLEMENTATION

  //! Virtual hook for create_ordinate_set
  void create_octant_ordinates_(std::vector<double> &mu, std::vector<double> &eta,
                                std::vector<double> &wt) const override;

  // DATA
};

} // end namespace rtt_quadrature

#endif // quadrature_Square_Chebyshev_Legendre_hh

//------------------------------------------------------------------------------------------------//
// end of quadrature/Square_Chebyshev_Legendre.hh
//------------------------------------------------------------------------------------------------//
