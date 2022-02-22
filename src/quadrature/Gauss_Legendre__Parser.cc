//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   quadrature/Gauss_Legendre.cc
 * \author Kelly Thompson
 * \date   Tue Feb 22 10:21:50 2000
 * \brief  A class representing an interval Gauss-Legendre quadrature set.
 * \note   Copyright (C) 2012-2022 Triad National Security, LLC., All rights reserved.
 *         reserved. */
//------------------------------------------------------------------------------------------------//

#include "Gauss_Legendre.hh"
#include "parser/utilities.hh"

namespace rtt_quadrature {
using namespace rtt_parser;

//------------------------------------------------------------------------------------------------//
std::shared_ptr<Quadrature> Gauss_Legendre::parse(Token_Stream &tokens) {
  Token token = tokens.shift();
  tokens.check_syntax(token.text() == "order", "expected an order");

  unsigned sn_order = parse_positive_integer(tokens);

  tokens.check_semantics(sn_order % 2 == 0, "order must be even");
  tokens.check_syntax(tokens.shift().type() == END, "missing end?");

  return std::make_shared<Gauss_Legendre>(sn_order);
}

} // end namespace rtt_quadrature

//------------------------------------------------------------------------------------------------//
// end of quadrature/Gauss_Legendre.cc
//------------------------------------------------------------------------------------------------//
