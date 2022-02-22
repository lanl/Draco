//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   quadrature/Gauss_Legendre.hh
 * \author Kelly Thompson
 * \date   Tue Feb 22 10:21:50 2000
 * \brief  A class representing an interval Gauss-Legendre quadrature set.
 * \note   Copyright (C) 2012-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_quadrature_Gauss_Legendre_hh
#define rtt_quadrature_Gauss_Legendre_hh

#include "Interval_Quadrature.hh"

namespace rtt_quadrature {

//================================================================================================//
/*!
 * \class Gauss_Legendre
 *
 * \brief A class representing an interval Gauss-Legendre quadrature set.
 *
 * This is an interval (e.g. 1D) angle quadrature set that achieves high formal accuracy by using
 * Gaussian integration based on the Legendre polynomials.
 */
//================================================================================================//
class Gauss_Legendre : public Interval_Quadrature {
public:
  // CREATORS
  explicit Gauss_Legendre(unsigned sn_order);

  // SERVICES

  std::string name() const override;
  std::string parse_name() const override;
  unsigned number_of_levels() const override;
  std::string as_text(std::string const &indent) const override;
  bool check_class_invariants() const;

  // STATICS

  static std::shared_ptr<Quadrature> parse(Token_Stream &tokens);

protected:
  std::vector<Ordinate> create_level_ordinates_(double norm) const override;
};

} // end namespace rtt_quadrature

#endif // rtt_quadrature_Quadrature_hh

//------------------------------------------------------------------------------------------------//
// end of quadrature/Quadrature.hh
//------------------------------------------------------------------------------------------------//
