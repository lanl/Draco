//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   parser/Text_Token_Stream.hh
 * \author Kent G. Budge
 * \brief  Definition of the Text_Token_Stream class.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_Text_Token_Stream_HH
#define rtt_Text_Token_Stream_HH

#include "Token_Stream.hh"
#include <set>
#include <stack>

namespace rtt_parser {
//------------------------------------------------------------------------------------------------//
/*!
 * \brief Abstract text-based token stream for simple parsers.
 *
 * A Text_Token_Stream constructs its stream of tokens by scanning a stream of text characters,
 * supplied through the \c getc method of its child classes.
 *
 * C and C++ comments are treated as breaking whitespace.
 *
 * Null characters are not permitted in the character stream.  They are used internally to indicate
 * the end of file or an error condition.
 */
class Text_Token_Stream : public Token_Stream {
public:
  // ACCESSORS

  //! Return the current line in the text stream.
  unsigned line() const {
    Ensure(line_ > 0);
    return line_;
  }

  bool no_nonbreaking_ws() const { return no_nonbreaking_ws_; }

  //! Return the current set of whitespace characters.
  std::set<char> const &whitespace() const { return whitespace_; }

  //! Check the class invariants.
  bool check_class_invariants() const;

  // MANIPULATORS

  void rewind() override = 0;

  // SERVICES

  //! Does the Token_Stream consider \c c to be whitespace?
  bool is_whitespace(char c) const;

  //! Does the Token_Stream consider \e c to be non-breaking
  bool is_nb_whitespace(char c) const;

  // CONST DATA

  //! The default whitespace definition
  DLL_PUBLIC_parser static std::set<char> const default_whitespace;

protected:
  // IMPLEMENTATION

  //! Construct a Text_Token_Stream.
  Text_Token_Stream();

  //! Construct a Text_Token_Stream.
  explicit Text_Token_Stream(std::set<char> const &ws, bool no_nonbreaking_ws = false);

  //! Scan the next token.
  Token fill_() override;

  //! Push a character onto the back of the character queue.
  void character_push_back_(char c);

  //! Move one or more characters from the text stream into the character buffer.
  virtual void fill_character_buffer_() = 0;

  virtual bool error_() const = 0;
  virtual bool end_() const = 0;
  virtual std::string location_() const = 0;

  //! Pop a character off the internal buffer.
  char pop_char_();
  //! Peek ahead at the internal buffer.
  char peek_(unsigned pos = 0);

  //! Skip any whitespace at the cursor position.
  void eat_whitespace_();

  //! Enter a nested file in a include directive.
  virtual void push_include(std::string &include_file_name) = 0;

  //! Exit a nested file from a include directive.
  virtual void pop_include() = 0;

  // The following scan_ functions are for numeric scanning.  The names reflect the context-free
  // grammar given by Stroustrup in appendix A of _The C++ Programming Language_.  However, we do
  // not presently recognize type suffixes on either integers or floats.
  unsigned scan_floating_literal_();
  unsigned scan_digit_sequence_(unsigned &pos);
  unsigned scan_exponent_part_(unsigned &pos);
  unsigned scan_fractional_constant_(unsigned &pos);

  unsigned scan_integer_literal_();
  unsigned scan_decimal_literal_(unsigned &pos);
  unsigned scan_hexadecimal_literal_(unsigned &pos);
  unsigned scan_octal_literal_(unsigned &pos);

  // Scan keyword
  Token scan_keyword();

  // Scan manifest string
  Token scan_manifest_string();

private:
  // IMPLEMENTATION

  // DATA

  std::stack<std::deque<char>> buffers_;
  std::deque<char> buffer_; //!< Character buffer. Refilled as needed using fill_character_buffer_()
  std::set<char> whitespace_;     //!< The whitespace character list
  std::stack<unsigned> lines_;    //!< Stack of current line values for nested input files.
  unsigned line_{1};              //!< Current line in input file.
  bool no_nonbreaking_ws_{false}; //!< Treat all whitespace as breaking whitespace.
};

} // namespace rtt_parser

#endif // rtt_Text_Token_Stream_HH
//------------------------------------------------------------------------------------------------//
// end of Text_Token_Stream.hh
//------------------------------------------------------------------------------------------------//
