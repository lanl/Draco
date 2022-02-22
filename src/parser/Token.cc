//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   Token.cc
 * \author Kent G. Budge
 * \brief  Definitions of Token helper functions.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

#include "Token.hh"
#include <cctype>
#include <cstdlib>
#include <cstring>

namespace rtt_parser {

//------------------------------------------------------------------------------------------------//
//! Is the argument a token type that has no associated text?
bool Is_Text_Token(Token_Type const type) {
  return type != rtt_parser::ERROR && type != EXIT && type != END;
}

//------------------------------------------------------------------------------------------------//
/*!
 * Is the argument a valid OTHER text?
 *
 * \return \c true if the argument points to a string of a single character that
 * does not fit any other token type pattern, or if the argument points to a
 * string of two or three characters from a recognized standard set.
 */
bool Is_Other_Text(char const *text) {
  Require(text != nullptr);

  if (text[0] == 0) {
    return false;
  } else if (text[1] == 0) {
    char const c = text[0];
    return !isalnum(c) && !isspace(c) && c != '_';
  } else if (text[2] == 0) {
    return !strcmp(text, "<=") || !strcmp(text, ">=") || !strcmp(text, "==") ||
           !strcmp(text, "!=") || !strcmp(text, "&&") || !strcmp(text, "||");
  } else
  // no three-character OTHER tokens recognized at present
  {
    return false;
  }
}

//------------------------------------------------------------------------------------------------//
/*!
 * Is the argument a valid keyword?
 *
 * \return \c true if the argument points to a string consisting of a sequence of C++ identifiers
 * separated by single spaces.
 */
bool Is_Keyword_Text(char const *text) {
  Require(text != nullptr);

  char c = *text++;
  while (true) {
    if (!isalpha(c) && c != '_')
      return false;
    while (c = *text++, isalnum(c) || c == '_') { /* do nothing */
    };
    if (!c)
      return true;
    if (c != ' ')
      return false;
    c = *text++;
  }
}

//------------------------------------------------------------------------------------------------//
/*!
 * Is the argument a valid string constant?
 *
 * \return \c true if the argument points to a string consisting of a single C++ string constant,
 * including the delimiting quotes.
 */
bool Is_String_Text(char const *text) {
  Require(text != nullptr);

  char c = *text++;
  if (c != '"')
    return false;
  while (true) {
    c = *text++;
    if (c == 0)
      return false;
    if (c == '"')
      return !*text++;
    if (c == '\\') {
      if (!*text++)
        return false;
    }
  }
}

//------------------------------------------------------------------------------------------------//
/*!
 *  Is the argument a valid real constant?
 *
 * \return \c true if the argument points to a string consisting of a single C++
 * floating-point constant.
 */
bool Is_Real_Text(char const *text) {
  Require(text != nullptr);

  char *endtext;
  strtod(text, &endtext);
  return endtext != text && *endtext == '\0';
}

//------------------------------------------------------------------------------------------------//
/*!
 * Is the argument a valid integer constanta?
 *
 * \return \c true if the argument points to a string consisting of a single C++
 * integer constant.
 */
bool Is_Integer_Text(char const *text) {
  Require(text != nullptr);

  char *endtext;
  strtol(text, &endtext, 0);
  return !*endtext;
}

//------------------------------------------------------------------------------------------------//
/*!
 * \param a First token to compare
 * \param b Second token to compare
 *
 * \return \c true if the two tokens are equal.
 */
bool operator==(Token const &a, Token const &b) {
  return a.type() == b.type() && a.text() == b.text() && a.location() == b.location();
}

//------------------------------------------------------------------------------------------------//
/*!
 * The invariants all reflect the basic requirement that the token text is
 * consistent with the token type.  For example, if the type is REAL, the text
 * must be a valid C representation of a real number, which can be converted to
 * double using atof.
 *
 * \return \c true if the invariants are all satisfied; \c false otherwise
 */
bool Token::check_class_invariant() const {
  return (Is_Text_Token(type_) || text_ == "") &&
         (type_ != KEYWORD || Is_Keyword_Text(text_.c_str())) &&
         (type_ != REAL || Is_Real_Text(text_.c_str())) &&
         (type_ != INTEGER || Is_Integer_Text(text_.c_str())) &&
         (type_ != STRING || Is_String_Text(text_.c_str())) &&
         (type_ != OTHER || Is_Other_Text(text_.c_str()));
}

} // namespace rtt_parser
//------------------------------------------------------------------------------------------------//
// end of Token_Stream.cc
//------------------------------------------------------------------------------------------------//
