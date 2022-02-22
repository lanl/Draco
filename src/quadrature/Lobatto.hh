//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   quadrature/Lobatto.hh
 * \author Kelly Thompson
 * \date   Tue Feb 22 10:21:50 2000
 * \brief  A class representing an interval Gauss-Legendre quadrature set.
 * \note   Copyright (C) 2012-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_quadrature_Lobatto_hh
#define rtt_quadrature_Lobatto_hh

#include "Interval_Quadrature.hh"

namespace rtt_quadrature {

//================================================================================================//
/*!
 * \class Lobatto
 *
 * \brief A class representing an interval Lobatto quadrature set.
 *
 * This is an interval (e.g. 1D) angle quadrature set whose abscissae at order N are the roots of
 * the derivative of the Laguerre polynomial of order N-1 plus the end points of the interval.
 */
//================================================================================================//

class Lobatto : public Interval_Quadrature {
public:
  // CREATORS
  explicit Lobatto(unsigned sn_order);

  // SERVICES

  std::string name() const override;
  std::string parse_name() const override;
  unsigned number_of_levels() const override;
  std::string as_text(std::string const &indent) const override;
  bool is_open_interval() const override;
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
