//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   parser/Token.hh
 * \author Kent G. Budge
 * \brief  Define class Token and enum Token_Type
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_Token_HH
#define rtt_Token_HH

#include "ds++/Assert.hh"

#if defined(MSVC)
#undef ERROR
#endif

namespace rtt_parser {
using std::string;

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Token types recognized by a Token_Stream.
 */
enum Token_Type {
  END,     /*!< The identifier <CODE>end</CODE>, denoting that the Parse_Table should return control
            *   to its client.  Can be used to implement nested parse tables.*/
  EXIT,    /*!< Denotes that the end of the Token_Stream has been reached.  The Token_Stream will
            *   continue to return EXIT indefinitely once its end has been reached. */
  KEYWORD, /*!< A sequence of one or more C++ identifiers separated by whitespace. */
  REAL,    /*!< A valid C++ floating-point constant. */
  INTEGER, /*!< A valid C++ integer constant. */
  STRING,  /*!< A valid C++ string constant. */
  ERROR,   /*!< The error token, indicating something wrong with the token stream.  For example, a
            *   file-based token stream would return this token if the file failed to open. */
  OTHER    /*!< A single character or sequence of characters (such as "==") that does not belong to
            *   one of the regular token types described above. */
};

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Description of a token.
 *
 * This class represents a lexical token, for use in simple parsing systems for analysis codes.  The
 * token is characterized by its type, value, and location.
 */
class Token {
public:
  // CREATORS

  //! Construct a Token with the specified non-text type and location.
  inline Token(Token_Type type, string const &location);

  //! Construct a single-character OTHER token with the specified location.
  inline Token(char c, string const &location);

  //! Construct a Token with specified type, text, and location.
  inline Token(Token_Type type, string const &text, string const &location);

  //! Construct a Token with specified type, text, and location.
  inline Token(Token_Type type, string &&text, string &&location);

  //! Default constructor
  inline Token() : text_(), location_() {}

  // ACCESSORS

  //! Return the token type.
  Token_Type type() const noexcept { return type_; }

  //! Return the token text.
  string const &text() const { return text_; }

  //! Return the location information.
  string const &location() const { return location_; }

  //! Check that the class invariants are satisfied.
  bool check_class_invariant() const;

  // MANIPULATORS

  void swap(Token &src) {
    std::swap(type_, src.type_);
    text_.swap(src.text_);
    location_.swap(src.location_);
  }

private:
  Token_Type type_ = {END}; //!< Type of this token
  string text_;             //!< Text of this token
  string location_;         //!< Location information (such as file and line)
};

// For checking of assertions
bool Is_Text_Token(Token_Type type);
bool Is_Integer_Text(char const *string);
bool Is_Keyword_Text(char const *string);
bool Is_Real_Text(char const *string);
bool Is_String_Text(char const *string);
bool Is_Other_Text(char const *string);

//! Test equality of two Tokens
bool operator==(Token const &a, Token const &b);

//------------------------------------------------------------------------------------------------//
/*!
 * \param type Type of the Token.
 * \param text Text of the Token.
 * \param location The token location.
 */
inline Token::Token(Token_Type const type, string const &text, string const &location)
    : type_(type), text_(text), location_(location) {
  Require(Is_Text_Token(type));
  Require(type != KEYWORD || Is_Keyword_Text(text.c_str()));
  Require(type != REAL || Is_Real_Text(text.c_str()));
  Require(type != INTEGER || Is_Integer_Text(text.c_str()));
  Require(type != STRING || Is_String_Text(text.c_str()));
  Require(type != OTHER || Is_Other_Text(text.c_str()));

  Ensure(check_class_invariant());
  Ensure(this->type() == type);
  Ensure(this->text() == text);
  Ensure(this->location() == location);
}

//------------------------------------------------------------------------------------------------//
/*!
 * Move version of previous constructor.
 *
 * \param type Type of the Token.
 * \param text Text of the Token.
 * \param location The token location.
 */
inline Token::Token(Token_Type const type, string &&text, string &&location)
    : type_(type), text_(text), location_(location) {
  Require(Is_Text_Token(type));
  Require(type != KEYWORD || Is_Keyword_Text(text.c_str()));
  Require(type != REAL || Is_Real_Text(text.c_str()));
  Require(type != INTEGER || Is_Integer_Text(text.c_str()));
  Require(type != STRING || Is_String_Text(text.c_str()));
  Require(type != OTHER || Is_Other_Text(text.c_str()));

  Ensure(check_class_invariant());
  Ensure(this->type() == type);
  Ensure(this->text() == text);
  Ensure(this->location() == location);
}

//------------------------------------------------------------------------------------------------//
/*!
 * \param c The token text (a single character)
 * \param location The token location.
 */
inline Token::Token(char const c, string const &location)
    : type_(OTHER), text_(1, c), location_(location) {
  Require(Is_Other_Text(string(1, c).c_str()));

  Ensure(check_class_invariant());
  Ensure(this->type() == OTHER);
  Ensure(this->text() == string(1, c));
  Ensure(this->location() == location);
}

//------------------------------------------------------------------------------------------------//
/*!
 * \param type Token type to create; must be one of END, EXIT, or ERROR.
 * \param location The token location
 *
 * These token types have no associated text.
 */
inline Token::Token(Token_Type const type, string const &location)
    : type_(type), text_(), location_(location) {
  Require(!Is_Text_Token(type));

  Ensure(check_class_invariant());
  Ensure(this->type() == type);
  Ensure(this->text() == "");
  Ensure(this->location() == location);
}

} // namespace rtt_parser

#endif // rtt_Token_HH

//------------------------------------------------------------------------------------------------//
// end of Token.hh
//------------------------------------------------------------------------------------------------//
