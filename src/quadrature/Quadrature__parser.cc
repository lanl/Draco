//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   quadrature/Quadrature__parser.cc
 * \author Kelly Thompson
 * \date   Tue Feb 22 10:21:50 2000
 * \brief  Parsers for various quadrature classes.
 * \note   Copyright (C) 2015-2022 Triad National Security, LLC., All rights reserved.
 *         reserved. */
//------------------------------------------------------------------------------------------------//

#include "Quadrature__parser.hh"

#include "Double_Gauss.hh"
#include "Gauss_Legendre.hh"
#include "General_Octant_Quadrature.hh"
#include "Level_Symmetric.hh"
#include "Lobatto.hh"
#include "Product_Chebyshev_Legendre.hh"
#include "Square_Chebyshev_Legendre.hh"
#include "Tri_Chebyshev_Legendre.hh"
#include "parser/Abstract_Class_Parser.hh"
#include "parser/utilities.hh"

using namespace rtt_parser;
using namespace rtt_quadrature;

namespace // anonymous
{
//------------------------------------------------------------------------------------------------//
std::shared_ptr<Quadrature> parse_gauss_legendre(Token_Stream &tokens) {
  return Gauss_Legendre::parse(tokens);
}

//------------------------------------------------------------------------------------------------//
std::shared_ptr<Quadrature> parse_level_symmetric(Token_Stream &tokens) {
  return Level_Symmetric::parse(tokens);
}

//------------------------------------------------------------------------------------------------//
std::shared_ptr<Quadrature> parse_tri_cl(Token_Stream &tokens) {
  return Tri_Chebyshev_Legendre::parse(tokens);
}

//------------------------------------------------------------------------------------------------//
std::shared_ptr<Quadrature> parse_square_cl(Token_Stream &tokens) {
  return Square_Chebyshev_Legendre::parse(tokens);
}

//------------------------------------------------------------------------------------------------//
std::shared_ptr<Quadrature> parse_product_cl(Token_Stream &tokens) {
  return Product_Chebyshev_Legendre::parse(tokens);
}

//------------------------------------------------------------------------------------------------//
std::shared_ptr<Quadrature> parse_double_gauss(Token_Stream &tokens) {
  return Double_Gauss::parse(tokens);
}

//------------------------------------------------------------------------------------------------//
std::shared_ptr<Quadrature> parse_general_octant_quadrature(Token_Stream &tokens) {
  return General_Octant_Quadrature::parse(tokens);
}

//------------------------------------------------------------------------------------------------//
std::shared_ptr<Quadrature> parse_lobatto(Token_Stream &tokens) { return Lobatto::parse(tokens); }

} // namespace

namespace rtt_parser {
//------------------------------------------------------------------------------------------------//
Class_Parse_Table<Quadrature> *Class_Parse_Table<Quadrature>::current_;
Parse_Table Class_Parse_Table<Quadrature>::parse_table_(nullptr, 0, Parse_Table::ONCE);
std::shared_ptr<Quadrature> Class_Parse_Table<Quadrature>::child_;

//------------------------------------------------------------------------------------------------//
Class_Parse_Table<Quadrature>::Class_Parse_Table() {
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
void Class_Parse_Table<Quadrature>::check_completeness(Token_Stream &tokens) {
  tokens.check_semantics(child_ != nullptr, "no quadrature specified");
}

//------------------------------------------------------------------------------------------------//
std::shared_ptr<Quadrature> Class_Parse_Table<Quadrature>::create_object() { return child_; }

//------------------------------------------------------------------------------------------------//
std::shared_ptr<Quadrature> &Class_Parse_Table<Quadrature>::get_parsed_object() { return child_; }

//------------------------------------------------------------------------------------------------//
/*!
 *
 * This function allows local developers to add new transport models to
 * Capsaicin, by adding their names to the set of transport models recognized by
 * the \c Transport_Model parser.
 *
 * \param keyword Name of the new transport model. Must be unique.
 * \param parse_function Function to call to parse a specification for the new
 * transport model. Must return an instance of the transport model class.
 */
/* static */
void Class_Parse_Table<Quadrature>::register_quadrature(
    string const &keyword, std::shared_ptr<Quadrature> parse_function(Token_Stream &)) {
  Abstract_Class_Parser<Quadrature, get_parse_table, get_parsed_object>::register_child(
      keyword, parse_function);
}

//------------------------------------------------------------------------------------------------//
template <> std::shared_ptr<Quadrature> parse_class<Quadrature>(Token_Stream &tokens) {
  Token token = tokens.shift();
  tokens.check_syntax(token.text() == "type", "expected type keyword");

  return parse_class_from_table<Class_Parse_Table<Quadrature>>(tokens);
}

} // namespace rtt_parser

//------------------------------------------------------------------------------------------------//
// end of quadrature/Quadrature.hh
//------------------------------------------------------------------------------------------------//
