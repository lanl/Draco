//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   parser/String_Token_Stream.hh
 * \author Kent G. Budge
 * \brief  Definition of class String_Token_Stream.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_String_Token_Stream_HH
#define rtt_String_Token_Stream_HH

#include "Text_Token_Stream.hh"
#include <fstream>

namespace rtt_parser {
using std::set;
using std::string;

//------------------------------------------------------------------------------------------------//
/*!
 * \brief std::string-based token stream
 *
 * String_Token_Stream is a Text_Token_Stream that obtains its text from a std::string passed to the
 * constructor. The diagnostic output is directed to an internal string that can be retrieved at
 * will.
 */
class String_Token_Stream : public Text_Token_Stream {
public:
  // CREATORS

  //! Construct a String_Token_Stream from a string (define a conversion path)
  String_Token_Stream(string const &text); // NOLINT [hicpp-explicit-conversions]

  //! Construct a String_Token_Stream from a string (define a conversion path)
  String_Token_Stream(string &&text); // NOLINT [hicpp-explicit-conversions]

  //! Construct a String_Token_Stream from a string.
  String_Token_Stream(string const &text, set<char> const &whitespace,
                      bool no_nonbreaking_ws = false);

  // MANIPULATORS

  // Return to the start of the string.
  void rewind() override;

  //! Report a condition.
  void report(Token const &token, string const &message) override;

  //! Report a condition.
  void report(string const &message) override;

  //! Report a comment.
  void comment(std::string const &message) override;

  // ACCESSORS

  //! Return the text to be tokenized.
  string const &text() const { return text_; }

  //! Return the accumulated set of messages.
  string messages() const { return messages_; }

  //! Check the class invariant.
  bool check_class_invariants() const;

protected:
  //! Generate a locator string.
  string location_() const override;

  void fill_character_buffer_() override;

  bool error_() const override;
  bool end_() const override;

  void push_include(std::string &include_file_name) override;
  void pop_include() override;

private:
  // IMPLEMENTATION

  // DATA

  string text_;  //!< Text to be tokenized
  unsigned pos_; //!< Cursor position in string

  string messages_; //!< Collection of diagnostic messages
};

} // namespace rtt_parser

#endif // rtt_String_Token_Stream_HH

//------------------------------------------------------------------------------------------------//
// end of String_Token_Stream.hh
//------------------------------------------------------------------------------------------------//
