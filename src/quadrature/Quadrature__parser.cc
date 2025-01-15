//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   quadrature/Quadrature__parser.cc
 * \author Kelly Thompson
 * \date   Tue Feb 22 10:21:50 2000
 * \brief  Parsers for various quadrature classes.
 * \note   Copyright (C) 2015-2025 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "quadrature/Quadrature__parser.hh"
#include "parser/Class_Parse_Table.hh"
#include "parser/Parse_Table.hh"
#include "parser/Token.hh"
#include "parser/Token_Stream.hh"
#include "parser/utilities.hh" // NOLINT(misc-include-cleaner)
#include "quadrature/Double_Gauss.hh"
#include "quadrature/Gauss_Legendre.hh"
#include "quadrature/General_Octant_Quadrature.hh"
#include "quadrature/Level_Symmetric.hh"
#include "quadrature/Lobatto.hh"
#include "quadrature/Product_Chebyshev_Legendre.hh"
#include "quadrature/Quadrature.hh"
#include "quadrature/Square_Chebyshev_Legendre.hh"
#include "quadrature/Tri_Chebyshev_Legendre.hh"
#include <memory>
#include <string>

//------------------------------------------------------------------------------------------------//
namespace // anonymous
{
//------------------------------------------------------------------------------------------------//
std::shared_ptr<rtt_quadrature::Quadrature> parse_gauss_legendre(rtt_parser::Token_Stream &tokens) {
  return rtt_quadrature::Gauss_Legendre::parse(tokens);
}

//------------------------------------------------------------------------------------------------//
std::shared_ptr<rtt_quadrature::Quadrature>
parse_level_symmetric(rtt_parser::Token_Stream &tokens) {
  return rtt_quadrature::Level_Symmetric::parse(tokens);
}

//------------------------------------------------------------------------------------------------//
std::shared_ptr<rtt_quadrature::Quadrature> parse_tri_cl(rtt_parser::Token_Stream &tokens) {
  return rtt_quadrature::Tri_Chebyshev_Legendre::parse(tokens);
}

//------------------------------------------------------------------------------------------------//
std::shared_ptr<rtt_quadrature::Quadrature> parse_square_cl(rtt_parser::Token_Stream &tokens) {
  return rtt_quadrature::Square_Chebyshev_Legendre::parse(tokens);
}

//------------------------------------------------------------------------------------------------//
std::shared_ptr<rtt_quadrature::Quadrature> parse_product_cl(rtt_parser::Token_Stream &tokens) {
  return rtt_quadrature::Product_Chebyshev_Legendre::parse(tokens);
}

//------------------------------------------------------------------------------------------------//
std::shared_ptr<rtt_quadrature::Quadrature> parse_double_gauss(rtt_parser::Token_Stream &tokens) {
  return rtt_quadrature::Double_Gauss::parse(tokens);
}

//------------------------------------------------------------------------------------------------//
std::shared_ptr<rtt_quadrature::Quadrature>
parse_general_octant_quadrature(rtt_parser::Token_Stream &tokens) {
  return rtt_quadrature::General_Octant_Quadrature::parse(tokens);
}

//------------------------------------------------------------------------------------------------//
std::shared_ptr<rtt_quadrature::Quadrature> parse_lobatto(rtt_parser::Token_Stream &tokens) {
  return rtt_quadrature::Lobatto::parse(tokens);
}

} // namespace

//------------------------------------------------------------------------------------------------//

namespace rtt_parser {

//------------------------------------------------------------------------------------------------//
Class_Parse_Table<rtt_quadrature::Quadrature>
    *Class_Parse_Table<rtt_quadrature::Quadrature>::current_;

Parse_Table Class_Parse_Table<rtt_quadrature::Quadrature>::parse_table_(nullptr, 0,
                                                                        Parse_Table::ONCE);

std::shared_ptr<rtt_quadrature::Quadrature> Class_Parse_Table<rtt_quadrature::Quadrature>::child_;

//------------------------------------------------------------------------------------------------//
Class_Parse_Table<rtt_quadrature::Quadrature>::Class_Parse_Table() {
  // Be sure parse table has standard models
  static bool first_time = true;
  if (first_time) {
    register_quadrature("gauss legendre", parse_gauss_legendre);
    register_quadrature("level symmetric", parse_level_symmetric);
    register_quadrature("tri cl", parse_tri_cl);
    register_quadrature("square cl", parse_square_cl);
    register_quadrature("product cl", parse_product_cl);
    register_quadrature("double gauss", parse_double_gauss);
    register_quadrature("general octant quadrature", parse_general_octant_quadrature);
    register_quadrature("lobatto", parse_lobatto);
    first_time = false;
  }
  current_ = this;
}

//------------------------------------------------------------------------------------------------//
void Class_Parse_Table<rtt_quadrature::Quadrature>::check_completeness(Token_Stream &tokens) {
  tokens.check_semantics(child_ != nullptr, "no quadrature specified");
}

//------------------------------------------------------------------------------------------------//
std::shared_ptr<rtt_quadrature::Quadrature>
Class_Parse_Table<rtt_quadrature::Quadrature>::create_object() {
  return child_;
}

//------------------------------------------------------------------------------------------------//
std::shared_ptr<rtt_quadrature::Quadrature> &
Class_Parse_Table<rtt_quadrature::Quadrature>::get_parsed_object() {
  return child_;
}

//------------------------------------------------------------------------------------------------//
template <>
std::shared_ptr<rtt_quadrature::Quadrature>
parse_class<rtt_quadrature::Quadrature>(Token_Stream &tokens) {
  Token const token = tokens.shift();
  tokens.check_syntax(token.text() == "type", "expected type keyword");

  return parse_class_from_table<Class_Parse_Table<rtt_quadrature::Quadrature>>(tokens);
}

} // namespace rtt_parser

//------------------------------------------------------------------------------------------------//
// end of quadrature/Quadrature.hh
//------------------------------------------------------------------------------------------------//
