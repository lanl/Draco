//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   parser/utilities.hh
 * \author Kent G. Budge
 * \brief  Declarations of a number of useful parsing utilities.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved.
 *
 * This file declares functions that parse certain common constructs in a uniform way.
 */
//------------------------------------------------------------------------------------------------//

#ifndef parser_utilities_hh
#define parser_utilities_hh

#include "Expression.hh"
#include "mesh_element/Geometry.hh"
#include <memory>

namespace rtt_parser {
//! Can the next token in the stream be interpreted as real number?
bool at_real(Token_Stream &tokens);

//! Is the next token in the stream a unit name?
bool at_unit_term(Token_Stream &tokens, unsigned position = 0);

unsigned parse_positive_integer(Token_Stream &tokens);

unsigned parse_unsigned_integer(Token_Stream &tokens);

int parse_integer(Token_Stream &tokens);

double parse_real(Token_Stream &tokens);

double parse_positive_real(Token_Stream &tokens);

double parse_nonnegative_real(Token_Stream &tokens);

bool parse_bool(Token_Stream &tokens);

Unit parse_unit(Token_Stream &tokens);

void parse_vector(Token_Stream &tokens, double *x);

//! parser a quote-delimited string, stripping the quotes.
std::string parse_manifest_string(Token_Stream &tokens);

void parse_geometry(Token_Stream &tokens, rtt_mesh_element::Geometry &parsed_geometry);

void parse_unsigned_vector(Token_Stream &tokens, unsigned *x, unsigned size);
void set_internal_unit_system(rtt_units::UnitSystem const &units);
void set_unit_expressions_are_required(bool b);
rtt_units::UnitSystem const &get_internal_unit_system();
bool unit_expressions_are_required();
void free_internal_unit_system();

//! parser a real number followed by a unit expression.
double parse_quantity(Token_Stream &tokens, Unit const &unit, char const *name);

//! parse an expression followed by a unit expression.
std::shared_ptr<Expression>
parse_quantity(Token_Stream &tokens, Unit const &unit, char const *name,
               unsigned number_of_variables,
               std::map<string, pair<unsigned, Unit>> const &variable_map);

double parse_temperature(Token_Stream &tokens);

std::shared_ptr<Expression>
parse_temperature(Token_Stream &tokens, unsigned number_of_variables,
                  std::map<string, pair<unsigned, Unit>> const &variable_map);

//------------------------------------------------------------------------------------------------//
/*! Template for parse function that produces a class object.
 *
 * The parse_class template function is intended for general use as a common signature for all
 * functions that parse an object of a specified class from a Token_Stream, so that the common code
 * idiom for parsing an object of MyClass is:
 *
 * \code
 *  auto spMyClass = parse_class<MyClass>(tokens);
 * \endcode
 *
 * Developers may specialize this function as needed. A particular implementation is suggested in
 * Class_Parse_Table.hh.
 *
 * \tparam class The type of the object for which a smart pointer is returned.
 * \tparam Context Parameter pack giving parsing context, if any.
 * \param tokens Token stream from which to parse the user input.
 * \param context Context for the parse.
 * \return A pointer to an object matching the user specification, or NULL if the specification is
 *         not valid.
 */
template <typename Class, typename... Context>
std::shared_ptr<Class> parse_class(Token_Stream &tokens, Context const &...context);

} // namespace rtt_parser

#endif
// parser_utilities_hh

//------------------------------------------------------------------------------------------------//
// end of utilities.hh
//------------------------------------------------------------------------------------------------//
