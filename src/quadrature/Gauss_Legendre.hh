//----------------------------------*-C++-*----------------------------------------------//
/*!
 * \file   quadrature/Gauss_Legendre.hh
 * \author Kelly Thompson
 * \date   Tue Feb 22 10:21:50 2000
 * \brief  A class representing an interval Gauss-Legendre quadrature set.
 * \note   Copyright 2016 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------------------//
// $Id: Quadrature.hh 6718 2012-08-30 20:03:01Z warsa $
//---------------------------------------------------------------------------------------//

#ifndef __quadrature_Gauss_Legendre_hh__
#define __quadrature_Gauss_Legendre_hh__

#include "Interval_Quadrature.hh"

namespace rtt_quadrature {

//=======================================================================================//
/*!
 * \class Gauss_Legendre
 *
 * \brief A class representing an interval Gauss-Legendre quadrature set.
 *
 * This is an interval (e.g. 1D) angle quadrature set that achieves high
 * formal accuracy by using Gaussian integration based on the Legendre polynomials.
 */
//=======================================================================================//

class Gauss_Legendre : public Interval_Quadrature {
public:
  // CREATORS
  DLL_PUBLIC_quadrature explicit Gauss_Legendre(unsigned sn_order);

  // ACCESSORS

  // SERVICES

  virtual string name() const;

  virtual string parse_name() const;

  virtual unsigned number_of_levels() const;

  virtual string as_text(string const &indent) const;

  bool check_class_invariants() const;

  // STATICS

  static SP<Quadrature> parse(Token_Stream &tokens);

protected:
  virtual vector<Ordinate> create_level_ordinates_(double norm) const;
};

} // end namespace rtt_quadrature

#endif // __quadrature_Quadrature_hh__

//---------------------------------------------------------------------------------------//
// end of quadrature/Quadrature.hh
//---------------------------------------------------------------------------------------//
