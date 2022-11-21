//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   String_Token_Stream.cc
 * \author Kent G. Budge
 * \date   Wed Jan 22 15:18:23 MST 2003
 * \brief  Definitions of String_Token_Stream methods.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "String_Token_Stream.hh"
#include "c4/C4_Functions.hh"
#include <sstream>

namespace rtt_parser {
using namespace std;
using namespace rtt_dsxx;

//------------------------------------------------------------------------------------------------//
/*!
 * Construct a String_Token_Stream that derives its text from the specified string. Use the default
 * Text_Token_Stream user-defined whitespace characters.
 *
 * \param text Text to be tokenized.
 */
String_Token_Stream::String_Token_Stream(string const &text)
    : text_(text), pos_(0), messages_(std::string()) {
  Ensure(check_class_invariants());
  Ensure(whitespace() == Text_Token_Stream::default_whitespace);
  Ensure(this->text() == text);
  Ensure(messages() == "");
}

//------------------------------------------------------------------------------------------------//
/*!
 * Construct a String_Token_Stream that derives its text from the specified string. Use the default
 * Text_Token_Stream user-defined whitespace characters.
 *
 * This is the move-aware version.
 *
 * \param text Text to be tokenized.
 */
String_Token_Stream::String_Token_Stream(string &&text)
    : text_(text), pos_(0), messages_(std::string()) {
  Ensure(check_class_invariants());
  Ensure(whitespace() == Text_Token_Stream::default_whitespace);
  Ensure(messages() == "");
}

//------------------------------------------------------------------------------------------------//
/*!
 * Construct a String_Token_Stream that derives its text from the specified string.
 *
 * \param text Text from which to extract tokens.
 * \param ws Points to a string containing user-defined whitespace characters.
 * \param no_nonbreaking_ws Causes spaces and tabs to be treated as breaking whitespace. This has
 *        the effect of forcing all keywords to consist of a single identifier.
 */
String_Token_Stream::String_Token_Stream(string const &text, set<char> const &ws,
                                         bool const no_nonbreaking_ws)
    : Text_Token_Stream(ws, no_nonbreaking_ws), text_(text), pos_(0), messages_(std::string()) {
  Ensure(check_class_invariants());
  Ensure(whitespace() == ws);
  Ensure(this->text() == text);
  Ensure(messages() == "");
}

//------------------------------------------------------------------------------------------------//
/*!
 * This function constructs and returns a string of the form "near \<text\>" where \<text\>
 * reproduces the region of text where the last token was parsed. This is useful for error reporting
 * in parsers.
 *
 * \return A string of the form "near <text>"
 */
string String_Token_Stream::location_() const {
  // search backwards four endlines
  unsigned begin;
  unsigned count = 0;
  for (begin = pos_; begin > 0; --begin) {
    if (text_[begin] == '\n') {
      if (++count == 4)
        break;
    }
  }
  Check(text_.size() < UINT_MAX);
  auto const end = static_cast<unsigned>(text_.size());
  unsigned i;
  for (i = begin; i < end; ++i) {
    char const c = text_[i];
    if (i >= pos_ && c == '\n') {
      break;
    }
  }
  // This kruftiness is to create the location string with a single allocation.
  string Result;
  Result.reserve(6 + i - begin);
  Result.insert(0U, "near\n", 5U);
  Result.insert(Result.end(), text_.begin() + begin, text_.begin() + i);
  Result.insert(Result.end(), 1U, '\n');
  return Result;
}

//------------------------------------------------------------------------------------------------//
void String_Token_Stream::fill_character_buffer_() {
  if (pos_ < text_.length()) {
    character_push_back_(text_[pos_++]);
  } else {
    character_push_back_('\x0');
  }

  Ensure(check_class_invariants());
}

//------------------------------------------------------------------------------------------------//
bool String_Token_Stream::error_() const { return false; }

//------------------------------------------------------------------------------------------------//
bool String_Token_Stream::end_() const { return pos_ >= text_.length(); }

//------------------------------------------------------------------------------------------------//
//! This function sends a messsage by writing it to an internal string.
void String_Token_Stream::report(Token const &token, string const &message) {
  messages_ += token.location() + "\n" + message + '\n';

  Ensure(check_class_invariants());
}

//------------------------------------------------------------------------------------------------//
/*!
 * This function sends a message by writing it to an internal string..
 *
 * This version assumes that the cursor is the error location.
 */
void String_Token_Stream::report(string const &message) {
  Token token = lookahead();
  messages_ += token.location() + "\n" + message + '\n';

  Ensure(check_class_invariants());
}

//------------------------------------------------------------------------------------------------//
/*!
 * This function sends a message by writing it to an internal string..
 *
 * This version prints no location information.
 */
void String_Token_Stream::comment(string const &message) {
  messages_ += message + '\n';

  Ensure(check_class_invariants());
}

//------------------------------------------------------------------------------------------------//
void String_Token_Stream::rewind() {
  pos_ = 0;

  Text_Token_Stream::rewind();

  Ensure(check_class_invariants());
}

//------------------------------------------------------------------------------------------------//
bool String_Token_Stream::check_class_invariants() const { return pos_ <= text_.length(); }

//------------------------------------------------------------------------------------------------//
/*!
 * String_Token_Stream does not presently support the include directive.
 */
void String_Token_Stream::push_include(std::string & /*include_file_name*/) {
  report_syntax_error("#include not supported for String_Token_Stream");
}

//------------------------------------------------------------------------------------------------//
/*!
 * String_Token_Stream does not presently support the include directive.
 */
void String_Token_Stream::pop_include() {
  /* this function should be unreachable. Please note this in code coverage. */
}

} // namespace rtt_parser

//------------------------------------------------------------------------------------------------//
// end of String_Token_Stream.cc
//------------------------------------------------------------------------------------------------//
