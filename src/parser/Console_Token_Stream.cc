//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   Console_Token_Stream.cc
 * \author Kent G. Budge
 * \brief  Definitions of Console_Token_Stream methods.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "Console_Token_Stream.hh"
#include <iostream>
#include <sstream>

namespace rtt_parser {
using namespace std;

//------------------------------------------------------------------------------------------------//
//! Use the default Text_Token_Stream user-defined whitespace characters.
Console_Token_Stream::Console_Token_Stream() {
  Ensure(check_class_invariants());
  Ensure(Console_Token_Stream::location_() == "input");
}

//------------------------------------------------------------------------------------------------//
/*!
 * \param ws User-defined whitespace characters.
 *
 * \param no_nonbreaking_ws Treat spaces and tabs as breaking whitespace. This has the effect of
 *           forcing all keywords to consist of a single identifier.
 */
Console_Token_Stream::Console_Token_Stream(set<char> const &ws, bool const no_nonbreaking_ws)
    : Text_Token_Stream(ws, no_nonbreaking_ws) {
  Ensure(check_class_invariants());
  Ensure(Console_Token_Stream::location_() == "input");
  Ensure(whitespace() == ws);
  Ensure(this->no_nonbreaking_ws() == no_nonbreaking_ws);
}

//------------------------------------------------------------------------------------------------//
/*!
 * For a Console_Token_Stream, location is not a terribly  meaningful concept. So we return "input"
 * as the location, which is true enough.
 *
 * \return The string "input".
 */
std::string Console_Token_Stream::location_() const { return "input"; }

//------------------------------------------------------------------------------------------------//
//! This function moves the next character from cin into the character buffer.
void Console_Token_Stream::fill_character_buffer_() {
  char c = static_cast<char>(cin.get());
  if (cin.fail()) {
    character_push_back_('\0');
  } else {
    if (c == '\n') {
      c = ';';
    }
    character_push_back_(c);
  }

  Ensure(check_class_invariants());
}

//------------------------------------------------------------------------------------------------//
bool Console_Token_Stream::error_() const { return cin.fail(); }

//------------------------------------------------------------------------------------------------//
/*!
 * This function may be used to check whether the user has typed an end of file character (ctrl-D on
 * most Unix systems).
 *
 * \return \c true if an end of file character has been typed; \c false otherwise.
 */

bool Console_Token_Stream::end_() const { return cin.eof(); }

//------------------------------------------------------------------------------------------------//
/*!
 * This function sends a message by writing it to the error console stream.
 */

void Console_Token_Stream::report(Token const &token, string const &message) {
  std::cerr << token.location() << ": " << message << std::endl;

  Ensure(check_class_invariants());
}

//------------------------------------------------------------------------------------------------//
/*!
 * This function sends a message by writing it to the error console stream. This version assumes
 * that the cursor gives the correct message location.
 */

void Console_Token_Stream::report(string const &message) {
  Token token = lookahead();
  std::cerr << token.location() << ": " << message << std::endl;

  Ensure(check_class_invariants());
}

//------------------------------------------------------------------------------------------------//
/*!
 * This function sends a message by writing it to the error console stream. This version prints no
 * location information.
 */

void Console_Token_Stream::comment(string const &message) {
  std::cerr << message << std::endl;

  Ensure(check_class_invariants());
}

//------------------------------------------------------------------------------------------------//
//! This function flushes cin and resets the error count.

void Console_Token_Stream::rewind() {
  cin.clear(); // Must clear the error/end flag bits.
  Text_Token_Stream::rewind();
  Ensure(check_class_invariants());
  Ensure(cin.rdstate() == 0);
  Ensure(location_() == "input");
}

//------------------------------------------------------------------------------------------------//
//! Console_Token_Stream does not presently support the include directive.
void Console_Token_Stream::push_include(std::string & /*include_file_name*/) {
  report_syntax_error("#include not supported for Console_Token_Stream");
}

//------------------------------------------------------------------------------------------------//
//! Console_Token_Stream does not presently support the include directive.
void Console_Token_Stream::pop_include() {
  /* this function should be unreachable. Please note this in code coverage. */
}

} // namespace rtt_parser

//------------------------------------------------------------------------------------------------//
// end of Console_Token_Stream.cc
//------------------------------------------------------------------------------------------------//
