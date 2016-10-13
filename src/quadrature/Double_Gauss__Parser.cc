//----------------------------------*-C++-*----------------------------------------------//
/*!
 * \file   quadrature/Double_Gauss.cc
 * \author Kelly Thompson
 * \date   Tue Feb 22 10:21:50 2000
 * \brief  A class representing an interval Gauss-Legendre quadrature set.
 * \note   Copyright 2016 Los Alamos National Security, LLC. All rights
 *         reserved. 
 */
//---------------------------------------------------------------------------------------//
// $Id: Quadrature.hh 6718 2012-08-30 20:03:01Z warsa $
//---------------------------------------------------------------------------------------//

#include "Double_Gauss.hh"
#include "parser/utilities.hh"

namespace rtt_quadrature {
using namespace rtt_parser;

//---------------------------------------------------------------------------------------//
/*static*/
SP<Quadrature> Double_Gauss::parse(Token_Stream &tokens) {
  Token token = tokens.shift();
  tokens.check_syntax(token.text() == "order", "expected an order");

  unsigned sn_order = parse_positive_integer(tokens);

  tokens.check_semantics(sn_order % 2 == 0, "order must be even");
  tokens.check_syntax(tokens.shift().type() == END, "missing end?");

  return SP<Quadrature>(new Double_Gauss(sn_order));
}

} // end namespace rtt_quadrature

//---------------------------------------------------------------------------------------//
//                       end of quadrature/Quadrature.cc
//---------------------------------------------------------------------------------------//
