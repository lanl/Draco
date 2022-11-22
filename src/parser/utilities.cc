//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   parser/utilities.cc
 * \author Kent G. Budge
 * \date   18 Feb 2003
 * \brief  Definitions of parsing utility functions.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "utilities.hh"
#include "units/PhysicalConstants.hh"
#include "units/PhysicalConstantsSI.hh"
#include <cctype>
#include <cerrno>
#include <cstdlib>
#include <sstream>

namespace // anonymous
{

bool require_unit_expressions = true;
// construct on first use (see
// https://isocpp.org/wiki/faq/ctors#static-init-order)
rtt_units::UnitSystem *internal_unit_system = nullptr;

} // namespace

namespace rtt_parser {
using namespace std;

//------------------------------------------------------------------------------------------------//
/*! Set the unit system to which all parser code converts unit expressions.
 *
 * By default, all unit expressions are converted to SI values by the parser utility routines. For
 * example, if the double parse_quantity function is given the text "2.7 cm", it will return a value
 * of 0.027. The function set_common_unit_system can be used to change the internal unit system, for
 * example, to CGS, so that if the double parse_quantity function is given the text "2.7 cm", it
 * will return a value of 2.7.
 *
 * Similarly, by default,the Expression parse_quantity function assumes that all input quantities to
 * the Expression it constructs will be in SI units, but this will be overridden if
 * set_common_unit_system is used to change the internal unit system. Thus, if Expression
 * parse_quantity is given the text "0.5*(t+2*x*s/cm)*erg" then by default it will construct an
 * Expression that computes 0.5e-7*(t+200*x). However, if set_internal_unit_system had previously
 * been called to set the internal unit system to CGS, the Expression constructed would compute
 * 0.5*(t+2*x).
 *
 * This function will normally be called once at the beginning of a parse, since if the internal
 * unit system changes halfway through a parse, there is rarely any easy way to go back and convert
 * quantities parsed earlier to the new internal unit system. For example, if an input file is
 * allowed to specify the internal unit system, this specification should normally be required to be
 * the first line in the input file.
 *
 * There is no delete associated with the new global static value. It remains in scope until the
 * program terminates.
 */
void set_internal_unit_system(rtt_units::UnitSystem const &units) {
  if (internal_unit_system != nullptr)
    delete internal_unit_system;
  internal_unit_system = new rtt_units::UnitSystem(units);
}

//---------------------------------------------------------------------------------------//
// For tracking down memory leaks, it is useful to be able to free any static data.

void free_internal_unit_system() {
  if (internal_unit_system != nullptr)
    delete internal_unit_system;
  internal_unit_system = nullptr;
}

//------------------------------------------------------------------------------------------------//
/*! Set whether unit expressions are mandatory.
 *
 * By default,the parse_quantity routines expect the texts they parse to be of the form "2.99792e10
 * cm/s", that is, a real number followed by a unit expression. However,
 * set_are_unit_expressions_required can be used to either turn on or turn off the mandatory unit
 * expression. Quantities can always have a unit expression, and the quantity will then be converted
 * to the internal unit system, but if the mandatory unit expression flag is turned off, quantities
 * may optionally omit the unit expression and will be assumed to be expressed in the internal unit
 * system.
 *
 * Thus, if unit expressions are mandatory and the internal unit system has been set to CGS, then
 * "2.7 m" will be read by double parse_quantity as 270 while "2.7" will be an error. If unit
 * expressions are not mandatory and the internal unit system has been set to CGS, then "2.7 m" will
 * still be read by double parse_quantity as 270, but "2.7" will be read as 2.7.
 *
 * It should be clear that mandatory unit expressions are preferable when humans are writing the
 * input files being parsed, since input files that consistently use unit expressions will be less
 * ambiguous, better documented, and more likely to be free from error. The ability to turn off
 * mandatory unit expressions is targeted primarily at situations where it is programs that are
 * communicating with each other through text files rather than humans communicating to a program.
 *
 * It is particularly not recommended that unit expressions be made optional where humans will be
 * composing Expressions, since the effect of turning off mandatory unit expressions is to turn off
 * checking of internal dimensional consistency of Expressions. This is a powerful error checking
 * capability that should not be discarded lightly.
 */
void set_unit_expressions_are_required(bool const b) { require_unit_expressions = b; }

//------------------------------------------------------------------------------------------------//
/*! Get the default unit system
 *
 * There is no delete associated with the new global static value. It remains in
 * scope until the program terminates.
 */
rtt_units::UnitSystem const &get_internal_unit_system() {
  if (internal_unit_system == nullptr)
    internal_unit_system = new rtt_units::UnitSystem();
  return *internal_unit_system;
}

//------------------------------------------------------------------------------------------------//
bool unit_expressions_are_required() { return require_unit_expressions; }

//------------------------------------------------------------------------------------------------//
/*!
 * \param[in,out] tokens Token stream from which to parse the quantity.
 *
 * \return The parsed quantity.
 */
unsigned parse_unsigned_integer(Token_Stream &tokens) {
  Token const token = tokens.shift();
  tokens.check_syntax(token.type() == INTEGER, "expected an unsigned integer");
  errno = 0;
  char *endptr;
  unsigned long const Result = strtoul(token.text().c_str(), &endptr, 0);

  tokens.check_semantics(Result == static_cast<unsigned>(Result) && errno != ERANGE,
                         "integer value overflows");

  Check(endptr != nullptr);
  // Check(Result<UINT32_MAX);
  return static_cast<unsigned>(Result);
}

//------------------------------------------------------------------------------------------------//
/*!
 * \param[in,out] tokens Token stream from which to parse the quantity.
 * \return The parsed quantity.
 */
unsigned parse_positive_integer(Token_Stream &tokens) {
  unsigned Result = parse_unsigned_integer(tokens);
  if (Result == 0) {
    tokens.report_semantic_error("expected a positive integer");
    Result = 1;
  }

  Ensure(Result > 0);
  return Result;
}

//------------------------------------------------------------------------------------------------//
/*!
 * \param[in,out] tokens Token stream from which to parse the quantity.
 * \return The parsed quantity.
 */
int parse_integer(Token_Stream &tokens) {
  Token token = tokens.shift();
  string text;
  if (token.text() == "+") {
    token = tokens.shift();
  } else if (token.text() == "-") {
    text = '-';
    token = tokens.shift();
  }
  if (token.type() == INTEGER) {
    text += token.text();
    errno = 0;
    char *endptr;
    const long Result = strtol(text.c_str(), &endptr, 0);
    if (Result != static_cast<int>(Result) || errno == ERANGE) {
      tokens.report_semantic_error("integer value overflows");
    }
    Check(endptr != nullptr && *endptr == '\0');
    // Check(std::abs(Result)<INT32_MAX);
    return static_cast<int>(Result);
  } else {
    tokens.report_syntax_error(token, "expected an integer");
    return 0;
  }
}

//------------------------------------------------------------------------------------------------//
/*!
 * This function does not move the cursor in the token stream.
 *
 * \param tokens Token stream from which to parse the quantity.
 * \return \c true if the next token is REAL or INTEGER; \c false otherwise.
 */
bool at_real(Token_Stream &tokens) {
  Token const &token = tokens.lookahead();
  if (token.text() == "-" || token.text() == "+") {
    Token const &token2 = tokens.lookahead(1);
    return (token2.type() == REAL || token2.type() == INTEGER);
  } else {
    return (token.type() == REAL || token.type() == INTEGER);
  }
}

//------------------------------------------------------------------------------------------------//
/*!
 * We permit an integer token to appear where a real is expected, consistent with the integers being
 * a subset of reals, and with about five decades of common practice in the computer language
 * community.
 *
 * \param[in,out] tokens Token stream from which to parse the quantity.
 * \return The parsed quantity.
 */
double parse_real(Token_Stream &tokens) {
  Token token = tokens.shift();
  double sign = 1.0;
  if (token.text() == "+") {
    token = tokens.shift();
  } else if (token.text() == "-") {
    sign = -1.0;
    token = tokens.shift();
  }
  if (token.type() == REAL || token.type() == INTEGER) {
    string const &text = token.text();
    errno = 0;
    char *endptr;
    const double Result = strtod(text.c_str(), &endptr);
    if (errno == ERANGE) {
      if (Result > 1.0 || Result < -1.0) {
        tokens.report_semantic_error("real value overflows: " + text);
      } else {
        tokens.report("note: real value underflows: " + text);
      }
    }
    Check(endptr != nullptr && *endptr == '\0');
    return sign * Result;
  } else {
    tokens.report_syntax_error(token, "expected a real number");
    return 0.0;
  }
}

//------------------------------------------------------------------------------------------------//
/*!
 * \param tokens Token stream from which to parse the quantity.
 * \return The parsed quantity.
 */
double parse_positive_real(Token_Stream &tokens) {
  double Result = parse_real(tokens);
  if (Result <= 0.0) {
    tokens.report_semantic_error("expected a positive quantity");
    Result = 1;
  }

  Ensure(Result > 0);
  return Result;
}

//------------------------------------------------------------------------------------------------//
/*!
 * \param tokens Token stream from which to parse the quantity.
 * \return The parsed quantity.
 */
double parse_nonnegative_real(Token_Stream &tokens) {
  double Result = parse_real(tokens);
  if (Result < 0.0) {
    tokens.report_semantic_error("expected a nonnegative quantity");
    Result = 1;
  }

  Ensure(Result >= 0);
  return Result;
}

//------------------------------------------------------------------------------------------------//
/*!
 * \param tokens Token stream from which to parse the quantity.
 * \param x On return, contains the parsed vector components.
 * \pre \c x!=nullptr
 */
void parse_vector(Token_Stream &tokens, double *x) {
  Require(x != nullptr);

  // At least one component must be present.
  x[0] = parse_real(tokens);

  if (at_real(tokens)) {
    x[1] = parse_real(tokens);
    if (at_real(tokens)) {
      x[2] = parse_real(tokens);
    } else {
      x[2] = 0.0;
    }
  } else {
    x[1] = 0.0;
    x[2] = 0.0;
  }
}

//------------------------------------------------------------------------------------------------//
/*!
 * \param tokens Token stream from which to parse the quantity.
 * \param x On return, contains the parsed vector components.
 * \param size size of parameter x.
 *
 * \pre \c x!=nullptr
 */
void parse_unsigned_vector(Token_Stream &tokens, unsigned *x, unsigned size) {
  Require(x != nullptr);

  for (unsigned i = 0; i < size; ++i) {
    if (at_real(tokens)) {
      x[i] = parse_unsigned_integer(tokens);
    } else {
      ostringstream str;
      str << "unsigned integer sequence too short; expected " << size << " values.";

      tokens.report_semantic_error(str.str().c_str());
      return;
    }
  }
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Are we at a unit term?
 *
 * We are at a unit term if the next token on the Token_Stream is a valid unit name.
 *
 * \param tokens Token_Stream from which to parse.
 * \param position Position in Token_Stream at which to parse.  This lookahead capability is needed
 *        by parse_unit to see if a hyphen '-' is part of a unit expression.
 *
 * \return \c true if we are at the start of a unit term; \c false otherwise
 */
bool at_unit_term(Token_Stream &tokens, unsigned position) {
  Token const &token = tokens.lookahead(position);
  if (token.type() == KEYWORD) {
    string const &u = token.text();
    switch (u[0]) {
    case 'A':
      // return (u[1]=='\0');
      return (u.size() == 1);

    case 'C':
      return (u.size() == 1);

    case 'F':
      return (u.size() == 1);

    case 'H':
      return (u.size() == 1 || token.text() == "Hz");

    case 'J':
      return (u.size() == 1);

    case 'K':
      return (u.size() == 1);

    case 'N':
      return (u.size() == 1);

    case 'P':
      return (token.text() == "Pa");

    case 'S':
      return (u.size() == 1);

    case 'T':
      return (u.size() == 1);

    case 'V':
      return (u.size() == 1);

    case 'W':
      return (u.size() == 1 || token.text() == "Wb");

    case 'c':
      return (token.text() == "cd" || token.text() == "cm");

    case 'd':
      return (token.text() == "dyne");

    case 'e':
      return (token.text() == "erg" || token.text() == "eV");

    case 'f':
      return (token.text() == "foot");

    case 'g':
      return (u.size() == 1);

    case 'i':
      return (token.text() == "inch");

    case 'j':
      return (token.text() == "jerk");

    case 'k':
      return (token.text() == "kg" || token.text() == "keV");

    case 'l':
      return (token.text() == "lm" || token.text() == "lx");

    case 'm':
      return (u.size() == 1 || token.text() == "mol");

    case 'o':
      return (token.text() == "ohm");

    case 'p':
      return (token.text() == "pound");

    case 'r':
      return (token.text() == "rad");

    case 's':
      return (u.size() == 1 || token.text() == "sr" || token.text() == "sh" ||
              token.text() == "shake");

    default:
      return false;
    }
  } else if (token.type() == OTHER && token.text() == "(") {
    return true;
  } else {
    return false;
  }
}

Unit parse_unit(Token_Stream &tokens);

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Parse a unit name.
 *
 * A unit name can either be a literal unit name, like "kg", or a parenthesized unit expression.
 *
 * \param tokens Token_Stream from which to parse.
 * \return The unit.
 */
static Unit parse_unit_name(Token_Stream &tokens) {
  // Return value
  Unit retval{};
  Token token = tokens.shift();
  if (token.type() == KEYWORD) {
    string const &u = token.text();
    switch (u[0]) {
    case 'A':
      if (u.size() == 1)
        retval = A;
      else
        tokens.report_syntax_error("expected a unit");
      break;

    case 'C':
      if (u.size() == 1)
        retval = C;
      else
        tokens.report_syntax_error("expected a unit");
      break;

    case 'F':
      if (u.size() == 1)
        retval = F;
      else
        tokens.report_syntax_error("expected a unit");
      break;

    case 'H':
      if (u.size() == 1)
        retval = H;
      else if (u.size() == 2)
        retval = Hz;
      else
        tokens.report_syntax_error("expected a unit");
      break;

    case 'J':
      if (u.size() == 1)
        retval = J;
      else
        tokens.report_syntax_error("expected a unit");
      break;

    case 'K':
      if (u.size() == 1)
        retval = K;
      else
        tokens.report_syntax_error("expected a unit");
      break;

    case 'N':
      if (u.size() == 1)
        retval = N;
      else
        tokens.report_syntax_error("expected a unit");
      break;

    case 'P':
      if (u.size() == 2)
        retval = Pa;
      else
        tokens.report_syntax_error("expected a unit");
      break;

    case 'S':
      if (u.size() == 1)
        retval = S;
      else
        tokens.report_syntax_error("expected a unit");
      break;

    case 'T':
      if (u.size() == 1)
        retval = T;
      else
        tokens.report_syntax_error("expected a unit");
      break;

    case 'V':
      if (u.size() == 1)
        retval = V;
      else
        tokens.report_syntax_error("expected a unit");
      break;

    case 'W':
      if (u.size() == 1)
        retval = W;
      else if (token.text() == "Wb")
        retval = Wb;
      else
        tokens.report_syntax_error("expected a unit");
      break;

    case 'c':
      if (token.text() == "cd")
        retval = cd;
      else if (token.text() == "cm")
        retval = cm;
      else
        tokens.report_syntax_error("expected a unit");
      break;

    case 'd':
      if (token.text() == "dyne")
        retval = dyne;
      else
        tokens.report_syntax_error("expected a unit");
      break;

    case 'e':
      if (token.text() == "erg")
        retval = erg;
      else if (token.text() == "eV")
        retval = eV;
      else
        tokens.report_syntax_error("expected a unit");
      break;

    case 'f':
      if (token.text() == "foot")
        retval = foot;
      else
        tokens.report_syntax_error("expected a unit");
      break;

    case 'g':
      if (u.size() == 1)
        retval = g;
      else
        tokens.report_syntax_error("expected a unit");
      break;

    case 'i':
      if (token.text() == "inch")
        retval = inch;
      else
        tokens.report_syntax_error("expected a unit");
      break;

    case 'j':
      if (token.text() == "jerk")
        retval = J * 1e9;
      else
        tokens.report_syntax_error("expected a unit");
      break;

    case 'k':
      if (token.text() == "kg")
        retval = kg;
      else if (token.text() == "keV")
        retval = keV;
      else
        tokens.report_syntax_error("expected a unit");
      break;

    case 'l':
      if (token.text() == "lm")
        retval = lm;
      else if (token.text() == "lx")
        retval = lx;
      else
        tokens.report_syntax_error("expected a unit");
      break;

    case 'm':
      if (u.size() == 1)
        retval = m;
      else if (token.text() == "mol")
        retval = mol;
      else
        tokens.report_syntax_error("expected a unit");
      break;

    case 'o':
      if (token.text() == "ohm")
        retval = ohm;
      else
        tokens.report_syntax_error("expected a unit");
      break;

    case 'p':
      if (token.text() == "pound")
        retval = pound;
      else
        tokens.report_syntax_error("expected a unit");
      break;

    case 'r':
      if (token.text() == "rad")
        retval = rad;
      else
        tokens.report_syntax_error("expected a unit");
      break;

    case 's':
      if (u.size() == 1)
        retval = s;
      else if (token.text() == "sr")
        retval = sr;
      else if (token.text() == "sh" || token.text() == "shake")
        retval = 1.0e-8 * s;
      else
        tokens.report_syntax_error("expected a unit");
      break;

    default:
      tokens.report_syntax_error("expected a unit");
    }
  } else if (token.type() == OTHER && token.text() == "(") {
    Unit Result = parse_unit(tokens);
    token = tokens.shift();
    if (token.type() != OTHER || token.text() != ")")
      tokens.report_syntax_error("missing ')'");
    retval = Result;
  } else {
    tokens.report_syntax_error("expected a unit expression");
  }
  // never reached but causes warnings
  // return dimensionless;
  return retval;
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Parse a unit term.
 *
 * A unit term is a unit name optionally raised to some power.
 *
 * \param tokens Token_Stream from which to parse.
 * \return The unit term.
 */
static Unit parse_unit_term(Token_Stream &tokens) {
  Unit const Result = parse_unit_name(tokens);
  Token const &token = tokens.lookahead();
  if (token.text() == "^") {
    tokens.shift();
    double const exponent = parse_real(tokens);
    return pow(Result, exponent);
  }
  return Result;
}

//------------------------------------------------------------------------------------------------//
/*!
 * A unit expression is a sequence of tokens with a form such as "kg-m/sec" or "erg/cm^2/sec/Hz"
 * that gives the dimensions of a physical quantity.  This function parses such an expression from
 * its Token_Stream, returning the result as a Unit whose conversion factor is relative to SI.  An
 * empty unit expression is permitted and returns rtt_parser::dimensionless, the identity Unit
 * representing the pure number 1.
 *
 * \param tokens Token_Stream from which to parse.
 * \return The unit expression.
 */
Unit parse_unit(Token_Stream &tokens) {
  if (!at_unit_term(tokens))
    return dimensionless;

  Unit Result = parse_unit_term(tokens);

  for (;;) {
    Token const &token = tokens.lookahead();
    if (token.type() == OTHER) {
      if (token.text() == "-" && at_unit_term(tokens, 1)) {
        tokens.shift();
        Result = Result * parse_unit_term(tokens);
      } else if (token.text() == "/") {
        tokens.shift();
        Result = Result / parse_unit_term(tokens);
      } else {
        return Result;
      }
    } else {
      return Result;
    }
  }
}

//------------------------------------------------------------------------------------------------//
/*!
 * This function parses a quantity having dimensions. It is assumed that the client expects certain
 * dimensions for the quantity, and an exception is thrown if the dimensions are not what the client
 * expected. The quantity will be converted to the desired unit system, as indicated by the \c .conv
 * member of the \c target_unit argument.
 *
 * \param tokens Token stream from which to parse the quantity.
 * \param target_unit Expected units for the quantity parsed, including* conversion factor.
 * \param name Name of the units expected for the quantity parsed, such as "length" or
 *        "ergs/cm/sec/Hz". Used to generate diagnostic messages.
 *
 * \return The parsed value, converted to the desired unit system.
 */
double parse_quantity(Token_Stream &tokens, Unit const &target_unit, char const *const name) {
  double const value = parse_real(tokens);
  if (!unit_expressions_are_required() && !at_unit_term(tokens)) {
    return value;
  } else {
    Unit const unit = parse_unit(tokens);
    if (!is_compatible(unit, target_unit)) {
      ostringstream buffer;
      buffer << "expected quantity with dimensions of " << name;
      tokens.report_semantic_error(buffer.str().c_str());
    }
    return value * unit.conv * conversion_factor(unit, get_internal_unit_system());
  }
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Parse a temperature specification
 *
 * It is very common for transport researchers to specify a temperature in units of energy, using
 * Boltzmann's constant as the conversion factor.  This function is useful for parsers that
 * accommodate this convention.
 *
 * \param tokens Token stream from which to parse the specification.
 * \return The parsed temperature.
 *
 * \post \c Result>=0.0
 */
double parse_temperature(Token_Stream &tokens) {
  double Temp = parse_real(tokens);

  if (Temp < 0.0) {
    tokens.report_semantic_error("temperature must be nonnegative");
    Temp = 0.0;
  }

  if (!unit_expressions_are_required() && !at_unit_term(tokens)) {
    return Temp;
  } else {
    Unit const u = parse_unit(tokens);
    Temp *= u.conv * conversion_factor(u, get_internal_unit_system());
    if (is_compatible(u, K)) {
      // no action needed
    } else if (is_compatible(u, J)) {
      Temp *= get_internal_unit_system().T() /
              (rtt_units::boltzmannSI * conversion_factor(u, get_internal_unit_system()));
    } else {
      tokens.report_semantic_error("expected quantity with units of "
                                   "temperature");
      return 0.0;
    }
    return Temp;
  }
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Parse a temperature specification
 *
 * It is very common for transport researchers to specify a temperature in units of energy, using
 * Boltzmann's constant as the conversion factor.  This function is useful for parsers that
 * accommodate this convention.
 *
 * This version of the function parses a temperature expression containing user-defined variables.
 *
 * \param tokens Token stream from which to parse the specification.
 * \param number_of_variables Number of user-defined variables potentially present in the parsed
 *        expression.
 * \param variable_map Map of variable names to variables indices.
 * \return The parsed temperature.
 *
 * \post \c Result>=0.0
 */
std::shared_ptr<Expression>
parse_temperature(Token_Stream &tokens, unsigned const number_of_variables,
                  std::map<string, pair<unsigned, Unit>> const &variable_map) {
  std::shared_ptr<Expression> Temp = Expression::parse(number_of_variables, variable_map, tokens);

  if (!unit_expressions_are_required() && !at_unit_term(tokens)) {
    return Temp;
  } else {
    rtt_units::UnitSystem const &unit_system = get_internal_unit_system();
    Unit const u = parse_unit(tokens) * Temp->units();
    double const conv = conversion_factor(u, unit_system);
    if (is_compatible(u, K)) {
      // no action needed
      Temp->set_units(conv * u);
    } else if (is_compatible(u, J)) {
      double const boltzmann = rtt_units::boltzmannSI * unit_system.e() / unit_system.T();

      Temp->set_units(conv * u / boltzmann);
    } else {
      tokens.report_syntax_error("expected quantity with units of "
                                 "temperature");
    }
    return Temp;
  }
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Parses a STRING token and strips the delimiting quotation marks.
 *
 * \param tokens Token_Stream from which to parse.
 * \return The stripped string.
 */
std::string parse_manifest_string(Token_Stream &tokens) {
  Token const token = tokens.shift();
  if (token.type() != STRING) {
    tokens.report_syntax_error("expected a string, but saw " + token.text());
  }
  string Result = token.text();
  string::size_type const length = Result.size();
  Check(length > 1);
  Result = Result.substr(1, length - 2);
  return Result;
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Parse a geometry specification.
 *
 * \param tokens Token stream from which to parse the geometry.
 * \param parsed_geometry On entry, if the value is not \c END_GEOMETRY, a diagnostic is generated
 *           to the token stream. On return, contains the geometry that was parsed.
 *
 * \post <code> parsed_geometry == rtt_mesh_element::AXISYMMETRIC || parsed_geometry ==
 *         rtt_mesh_element::CARTESIAN || parsed_geometry == rtt_mesh_element::SPHERICAL </code>
 */
void parse_geometry(Token_Stream &tokens, rtt_mesh_element::Geometry &parsed_geometry) {
  tokens.check_semantics(parsed_geometry == rtt_mesh_element::Geometry::END_GEOMETRY,
                         "geometry specified twice");

  Token const token = tokens.shift();
  if (token.text() == "axisymmetric" || token.text() == "cylindrical") {
    parsed_geometry = rtt_mesh_element::Geometry::AXISYMMETRIC;
  } else if (token.text() == "cartesian" || token.text() == "xy" || token.text() == "slab") {
    parsed_geometry = rtt_mesh_element::Geometry::CARTESIAN;
  } else if (token.text() == "spherical") {
    parsed_geometry = rtt_mesh_element::Geometry::SPHERICAL;
  } else {
    tokens.report_syntax_error(token, "expected a geometry option, but saw " + token.text());
  }
  Ensure(parsed_geometry == rtt_mesh_element::Geometry::AXISYMMETRIC ||
         parsed_geometry == rtt_mesh_element::Geometry::CARTESIAN ||
         parsed_geometry == rtt_mesh_element::Geometry::SPHERICAL);
  return;
}

//------------------------------------------------------------------------------------------------//
/*!
 * \param tokens Token stream from which to parse the quantity.
 * \return The parsed quantity.
 */
bool parse_bool(Token_Stream &tokens) {
  Token const token = tokens.shift();
  if (token.text() == "true") {
    return true;
  } else {
    tokens.check_syntax(token.text() == "false", "expected 'true' or 'false'");
    return false;
  }
}

//------------------------------------------------------------------------------------------------//
/*!
 * This function parses an expression having dimensions. It is assumed that the client expects
 * certain dimensions for the quantity, and an exception is thrown if the dimensions are not what
 * the client expected. The quantity will be converted to the desired unit system, as indicated by
 * the \c .conv member of the \c target_unit argument.
 *
 * \param tokens Token stream from which to parse the quantity.
 * \param target_unit Expected units for the quantity parsed, including conversion factor.
 * \param name Name of the units expected for the quantity parsed, such as "length" or
 *        "ergs/cm/sec/Hz". Used to generate diagnostic messages.
 * \param number_of_variables Number of user-defined variables potentially present in the parsed
 &        expression.
 * \param variable_map Map of variable names to variables indices.
 * \return The parsed value, converted to the desired unit system.
 */
std::shared_ptr<Expression>
parse_quantity(Token_Stream &tokens, Unit const &target_unit, char const *const name,
               unsigned const number_of_variables,
               std::map<string, pair<unsigned, Unit>> const &variable_map) {
  std::shared_ptr<Expression> value = Expression::parse(number_of_variables, variable_map, tokens);

  if (!unit_expressions_are_required() && !at_unit_term(tokens)) {
    return value;
  } else {
    Unit unit = parse_unit(tokens) * value->units();
    if (!is_compatible(unit, target_unit)) {
      ostringstream buffer;
      buffer << "expected quantity with dimensions of " << name;
      tokens.report_semantic_error(buffer.str().c_str());
    }
    unit.conv *= conversion_factor(unit, get_internal_unit_system());
    value->set_units(unit);
    return value;
  }
}

} // namespace rtt_parser

//------------------------------------------------------------------------------------------------//
// end of utilities.cc
//------------------------------------------------------------------------------------------------//
