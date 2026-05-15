//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   quadrature/Level_Symmetric.cc
 * \author Kelly Thompson
 * \brief  A class representing an interval Gauss-Legendre quadrature set.
 * \note   Copyright (C) 2012-2026 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "parser/Token.hh"
#include "parser/Token_Stream.hh"
#include "parser/utilities.hh"
#include "quadrature/Level_Symmetric.hh"
#include "quadrature/Quadrature.hh"
#include <memory> // shared_ptr

namespace rtt_quadrature {
using namespace rtt_parser;

//------------------------------------------------------------------------------------------------//
std::shared_ptr<Quadrature> Level_Symmetric::parse(rtt_parser::Token_Stream &tokens) {
  Token const token = tokens.shift();
  tokens.check_syntax(token.text() == "order", "expected an order");

  unsigned const sn_order = parse_positive_integer(tokens);
  tokens.check_semantics(sn_order % 2 == 0, "order must be even");

  tokens.check_syntax(tokens.shift().type() == END, "missing end?");

  return std::shared_ptr<Quadrature>(new Level_Symmetric(sn_order));
}

} // end namespace rtt_quadrature

//------------------------------------------------------------------------------------------------//
// end of quadrature/Quadrature.cc
//------------------------------------------------------------------------------------------------//
