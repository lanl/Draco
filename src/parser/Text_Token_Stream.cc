//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   parser/Text_Token_Stream.cc
 * \author Kent G. Budge
 * \brief  Contains definitions of all Text_Token_Stream member functions.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "Text_Token_Stream.hh"
#include "ds++/path.hh"
#include <cctype>
#include <cstring>
#include <string>

#if defined(MSVC)
#undef ERROR
#endif

namespace rtt_parser {
using namespace std;

//------------------------------------------------------------------------------------------------//
// Helper function to allow a string && argument to take over storage of a string.
static string give(string &source) {
  string Result;
  Result.swap(source);
  return Result;
}

//------------------------------------------------------------------------------------------------//
set<char> const Text_Token_Stream::default_whitespace = {'=', ':', ';', ','};

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Constructs a Text_Token_Stream with the specified set of breaking whitespace characters.
 *
 * \param ws String containing the user-defined whitespace characters for this Text_Token_Stream.
 *
 * \param no_nonbreaking_ws If true, treat spaces and tabs as breaking whitespace. This has the
 * effect of forcing all keywords to consist of a single identifier.
 *
 * Whitespace characters are classified as breaking or nonbreaking whitespace.  Nonbreaking
 * whitespace separates non-keyword tokens and identifiers within a keyword but has no other
 * significance. Breaking whitespace is similar to nonbreaking whitespace except that it always
 * separates tokens; thus, two identifiers separated by breaking whitespace are considered to belong
 * to separate keywords.
 *
 * Nonbreaking whitespace characters are the space and horizontal tab characters.
 *
 * Breaking whitespace characters include all other characters for which the standard C library
 * function <CODE>isspace(char)</CODE> returns a nonzero value, plus additional characters defined
 * as nonbreaking whitespace by the client of the Token_Stream. In particular, a newline character
 * is always breaking whitespace.
 *
 * Whitespace is stripped from the beginning and end of every token, and the nonbreaking whitespace
 * separating each identifier within a keyword is replaced by a single space character.
 */
Text_Token_Stream::Text_Token_Stream(set<char> const &ws, bool const no_nonbreaking_ws)
    : buffer_(), whitespace_(ws), no_nonbreaking_ws_(no_nonbreaking_ws) {
  Ensure(check_class_invariants());
  Ensure(ws == whitespace());
  Ensure(line() == 1);
  Ensure(this->no_nonbreaking_ws() == no_nonbreaking_ws);
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Constructs a Text_Token_Stream with the default set of breaking whitespace characters. See
 * the previous constructor documentation for a discussion of how whitespace is defined and handled.
 *
 * The default whitespace characters are contained in the set
 * \c Text_Token_Stream::default_whitespace.
 */
Text_Token_Stream::Text_Token_Stream() : buffer_(), whitespace_(default_whitespace) {
  Ensure(check_class_invariants());
  Ensure(whitespace() == default_whitespace);
  Ensure(line() == 1);
  Ensure(!no_nonbreaking_ws());
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Scan the next token from the character stream. The character stream is accessed via the
 * fill_character_buffer, error, and end functions, which are pure virtual functions.
 */
Token Text_Token_Stream::fill_() {
  eat_whitespace_();

  char c = peek_(); // Character at the current cursor position

  string token_location = location_();

  if (c == '\0') {
    // Sentinel value for error or end of file.
    if (end_()) {
      if (lines_.size() > 0) {
        // We're at the end of an included file. Pop up.
        pop_include();
        return fill_();
      } else {
        // We're at the end of the top file.
        Ensure(check_class_invariants());
        return {EXIT, token_location};
      }
    } else {
      // Code only reachable if hardware failure occurs. Code coverage not expected.
      Ensure(check_class_invariants());
      return {rtt_parser::ERROR, token_location};
    }
  } else {
    if (isalpha(c) || c == '_')
    // Beginning of a keyword or END token
    {
      Token Result = scan_keyword();
      Ensure(check_class_invariants());
      return Result;
    } else if (isdigit(c) || c == '.') {
      // A number of some kind.  Note that an initial sign ('+' or '-') is tokenized independently,
      // because it could be interpreted as a binary operator in arithmetic expressions.  It is up
      // to the parser to decide if this is the correct interpretation.
      unsigned const float_length = scan_floating_literal_();
      unsigned const int_length = scan_integer_literal_();
      string text;
      if (float_length > int_length) {
        text.reserve(float_length);
        for (unsigned i = 0; i < float_length; i++) {
          c = pop_char_();
          text += c;
        }
        Ensure(check_class_invariants());
        return {REAL, give(text), give(token_location)};
      } else if (int_length > 0) {
        text.reserve(int_length);
        for (unsigned i = 0; i < int_length; i++) {
          char cc = pop_char_();
          text += cc;
        }
        Ensure(check_class_invariants());
        return {INTEGER, give(text), give(token_location)};
      } else {
        Check(c == '.');
        pop_char_();
        Ensure(check_class_invariants());
        return {'.', give(token_location)};
      }
    } else if (c == '"')
    // Manifest string
    {
      Token Result = scan_manifest_string();
      Ensure(check_class_invariants());
      return Result;
    } else if (c == '#')
    // #directive
    {
      pop_char_();
      eat_whitespace_();
      c = peek_();
      if (!isalpha(c) && c != '_') {
        report_syntax_error("ill-formed #directive");
      } else {
        Token directive = scan_keyword();
        if (directive.text() == "include") {
          eat_whitespace_();
          if (peek_() == '"') {
            Token file = scan_manifest_string();
            string name = file.text();
            // strip end quotes. May allow internal quotes someday ...
            name = name.substr(1, name.size() - 2);
            push_include(name);
            return fill_();
          } else {
            report_syntax_error("#include requires file name in quotes");
          }
        } else {
          report_syntax_error("unrecognized #directive: #" + directive.text());
        }
      }
    } else if (c == '<')
    // Multicharacter OTHER
    {
      pop_char_();
      if (peek_() == '=') {
        pop_char_();
        Ensure(check_class_invariants());
        return {OTHER, "<=", give(token_location)};
      } else {
        Ensure(check_class_invariants());
        return {c, give(token_location)};
      }
    } else if (c == '>')
    // Multicharacter OTHER
    {
      pop_char_();
      if (peek_() == '=') {
        pop_char_();
        Ensure(check_class_invariants());
        return {OTHER, ">=", give(token_location)};
      } else {
        Ensure(check_class_invariants());
        return {c, give(token_location)};
      }
    } else if (c == '&')
    // Multicharacter OTHER
    {
      pop_char_();
      if (peek_() == '&') {
        pop_char_();
        Ensure(check_class_invariants());
        return {OTHER, "&&", give(token_location)};
      } else {
        Ensure(check_class_invariants());
        return {c, give(token_location)};
      }
    } else if (c == '|')
    // Multicharacter OTHER
    {
      pop_char_();
      if (peek_() == '|') {
        pop_char_();
        Ensure(check_class_invariants());
        return {OTHER, "||", give(token_location)};
      } else {
        Ensure(check_class_invariants());
        return {c, give(token_location)};
      }
    } else {
      // OTHER
      pop_char_();
      Ensure(check_class_invariants());
      return {c, give(token_location)};
    }
  }
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief This function searches for the argument character in its internal list of whitespace
 *        characters.
 *
 * \param c Character to be checked against the whitespace list.
 *
 * \return \c true if and only if the character is found in the internal whitespace list.
 */
bool Text_Token_Stream::is_whitespace(char const c) const {
  return isspace(c) || whitespace_.count(c);
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief This function searches for the argument character in the Token_Stream's internal list of
 * nonbreaking whitespace characters.
 *
 * \param c Character to be checked against the nonbreaking whitespace list.
 *
 * \return \c true if and only if the character is found in the internal nonbreaking whitespace
 * list, and is \e not found in the breaking whitespace list..
 */
bool Text_Token_Stream::is_nb_whitespace(char const c) const {
  return !whitespace_.count(c) && (c == ' ' || c == '\t');
}

//------------------------------------------------------------------------------------------------//
/*!
 * An internal buffer is used to implement unlimited lookahead, necessary for scanning numbers
 * (which have a quite complex regular expression.)  This function pops a character off the top of
 * the internal buffer, using fill_character_buffer() if necessary to ensure that there is at least
 * one character in the buffer.  If the next character is a carriage return, the line count is
 * incremented.
 *
 * \return The next character in the buffer.
 */
char Text_Token_Stream::pop_char_() {
  Remember(unsigned const old_line = line_);

  char const Result = peek_();
  buffer_.pop_front();
  if (Result == '\n') {
    line_++;
  }

  Ensure(check_class_invariants());
  Ensure((Result == '\n' && line_ == old_line + 1) || line_ == old_line);
  return Result;
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Try to scan a floating literal.
 *
 * \return length of scanned literal; 0 if no literal could be scanned.
 */
unsigned Text_Token_Stream::scan_floating_literal_() {
  unsigned pos = 0;
  if (scan_fractional_constant_(pos) > 0) {
    scan_exponent_part_(pos);
    return pos;
  } else if (scan_digit_sequence_(pos)) {
    if (scan_exponent_part_(pos) == 0)
      return 0;
    return pos;
  } else {
    return 0;
  }
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Try to scan a digit sequence.
 *
 * \return length of scanned text; 0 if no text could be scanned.
 */
unsigned Text_Token_Stream::scan_digit_sequence_(unsigned &pos) {
  unsigned const old_pos = pos;
  while (isdigit(peek_(pos)))
    pos++;
  return pos - old_pos;
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Try to scan an exponent part.
 *
 * \return length of scanned text; 0 if no text could be scanned.
 */
unsigned Text_Token_Stream::scan_exponent_part_(unsigned &pos) {
  unsigned const old_pos = pos;
  char c = peek_(pos);
  if (c == 'e' || c == 'E') {
    pos++;
    c = peek_(pos);
    if (c == '-' || c == '+')
      pos++;
    if (!scan_digit_sequence_(pos)) {
      pos = old_pos;
    }
  }
  return pos - old_pos;
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Try to scan a fractional constant.
 *
 * \return length of scanned text; 0 if no text could be scanned.
 */
unsigned Text_Token_Stream::scan_fractional_constant_(unsigned &pos) {
  unsigned const old_pos = pos;
  if (scan_digit_sequence_(pos) > 0) {
    if (peek_(pos) != '.') {
      pos = old_pos;
    } else {
      pos++;
      scan_digit_sequence_(pos);
    }
  } else if (peek_(pos) == '.') {
    pos++;
    if (scan_digit_sequence_(pos) == 0) {
      pos = old_pos;
    }
  }
  return pos - old_pos;
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Try to scan an integer literal.
 *
 * \return length of scanned literal; 0 if no literal could be scanned.
 */
unsigned Text_Token_Stream::scan_integer_literal_() {
  unsigned pos = 0;
  if (scan_decimal_literal_(pos) > 0) {
  } else if (scan_hexadecimal_literal_(pos)) {
  } else if (scan_octal_literal_(pos)) {
  } else {
    return 0;
  }
  return pos;
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Try to scan decimal literal.
 *
 * \return length of scanned text; 0 if no text could be scanned.
 */
unsigned Text_Token_Stream::scan_decimal_literal_(unsigned &pos) {
  unsigned const old_pos = pos;
  char c = peek_(pos);
  if (isdigit(c) && c != '0') {
    while (isdigit(c)) {
      pos++;
      c = peek_(pos);
    }
  }
  return pos - old_pos;
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Try to scan hexadecimal literal.
 *
 * \return length of scanned text; 0 if no text could be scanned.
 */
unsigned Text_Token_Stream::scan_hexadecimal_literal_(unsigned &pos) {
  unsigned old_pos = pos;
  if (peek_(pos) == '0') {
    pos++;
    char c = peek_(pos);
    if (c == 'x' || c == 'X') {
      pos++;
      c = peek_(pos);
      if (!isxdigit(c)) {
        pos = old_pos;
      } else {
        while (isxdigit(c)) {
          pos++;
          c = peek_(pos);
        }
      }
    } else {
      pos = old_pos;
    }
  }
  return pos - old_pos;
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Try to scan octal literal.
 *
 * \return length of scanned text; 0 if no text could be scanned.
 */
unsigned Text_Token_Stream::scan_octal_literal_(unsigned &pos) {
  unsigned old_pos = pos;
  char c = peek_(pos);
  if (c == '0') {
    while (isdigit(c) && c != '8' && c != '9') {
      pos++;
      c = peek_(pos);
    }
  }
  return pos - old_pos;
}

//------------------------------------------------------------------------------------------------//
/*!
 * An internal buffer is used to implement unlimited lookahead, necessary for scanning numbers
 * (which have a quite complex regular expression.)  This function peeks ahead the specified number
 * of places in the buffer, using fill_character_buffer() if necessary to ensure that there is a
 * sufficient number of characters in the buffer.
 *
 * \param pos Position at which to peek.
 *
 * \return The character at \c buffer[pos].
 */
char Text_Token_Stream::peek_(unsigned const pos) {
  while (buffer_.size() <= pos) {
    fill_character_buffer_();
  }

  Ensure(check_class_invariants());
  return buffer_[pos];
}

//------------------------------------------------------------------------------------------------//
/*!
 * This function flushes the Text_Token_Stream's internal buffers, so that scanning resumes at the
 * beginning of the file stream.  It is normally called by its overriding version in children of
 * Text_Token_Stream.
 */
void Text_Token_Stream::rewind() {
  while (!buffers_.empty())
    buffers_.pop();
  buffer_.clear();
  while (!lines_.empty())
    lines_.pop();
  line_ = 1;

  Token_Stream::rewind();

  Ensure(check_class_invariants());
  Ensure(error_count() == 0);
}

//------------------------------------------------------------------------------------------------//
bool Text_Token_Stream::check_class_invariants() const { return line_ > 0; }

//------------------------------------------------------------------------------------------------//
/*!
 * This function skips past any whitespace present at the cursor position, leaving the cursor at the
 * first non-whitespace character following the initial cursor position.
 */

/* private */
void Text_Token_Stream::eat_whitespace_() {
  for (;;) {
    // Scan whitespace
    char c = peek_();
    while (is_whitespace(c) && c != '\0') {
      pop_char_();
      c = peek_();
    }

    // Check for a comment
    if (c == '/') {
      if (peek_(1) == '/') {
        // C++ comment
        while (c != '\n' && !error_() && !end_()) {
          pop_char_();
          c = peek_();
        }
      } else if (peek_(1) == '*') {
        pop_char_(); // pop the '/'
        pop_char_(); // pop the '*'
        while ((peek_(0) != '*' || peek_(1) != '/') && !error_() && !end_()) {
          pop_char_();
        }
        pop_char_(); // pop the '*'
        pop_char_(); // pop the '/'
      } else {
        break;
      }
    } else {
      break;
    }
  }
  // private member function -- no invariant check
}

//------------------------------------------------------------------------------------------------//
/*!
 * \param c Character to be pushed onto the back of the character queue.
 */
void Text_Token_Stream::character_push_back_(char const c) {
  Remember(auto const old_buffer_size = static_cast<unsigned>(buffer_.size()));
  buffer_.push_back(c);
  Ensure(check_class_invariants());
  Ensure(buffer_.size() == old_buffer_size + 1);
  Ensure(buffer_.back() == c);
}

//------------------------------------------------------------------------------------------------//
/*!
 * \param file_name Name of file to be included at this point. On return, replaced with an absolute
 * path based on DRACO_INCLUDE_PATH if the relative path did not exist.
 *
 * Child classes must set a policy on how to include a file. For example, a File_Token_Stream can be
 * expected to read the included file in serial; a Parallel_File_Token_Stream can be expected to
 * read the included file in parallel; and a Console_Token_Stream or String_Token_Stream presently
 * do not provide this capability and will treat a include directive as an error.
 *
 * This function is pure virtual with an implementation. This means that every child class must
 * implement this function, but part of its implementation must be to include
 *
 * \code
 *    Text_Token_Stream::push_include(include_file_name);
 * \endcode
 *
 * as the first line in its implementation of this function. This call stashes the line and
 * character buffer of the underlying Text_Token_Stream and also finds the absolute path of the file
 * name.
 */
void Text_Token_Stream::push_include(std::string &file_name) {
  lines_.push(line_);
  line_ = 1;
  buffers_.push(buffer_);
  buffer_.clear();

  // Now find the absolute path of the file name.

  if (!rtt_dsxx::fileExists(file_name)) {
    // File name as given does not resolve to an existing file. Assume this is a relative path and
    // look for the file in DRACO_INCLUDE_PATH.

    // At present, DRACO_INCLUDE_PATH can contain only a single path.  Multiple search options may
    // be implemented in the future.
    auto path = getenv("DRACO_INCLUDE_PATH");
    if (path != nullptr) {
      // For now, the only other possibility:
      file_name = path + ('/' + file_name);
      // If this doesn't exist, the stream will report the error later on.
    }
    // else return and let the stream report the error
  }
  // else the file name as given resolves to an existing file

  Require(check_class_invariants());
  Require(line() == 1);
}

//------------------------------------------------------------------------------------------------//
/*!
 * This function is pure virtual with an implementation. This means that every child class must
 * implement this function, but part of its implementation must be to reset the line number by
 * directly calling the base version. That is, every child class must include
 *
 * \code
 *    Text_Token_Stream::pop_include(include_file_name);
 * \endcode
 *
 * as the first line in its implementation of this function.
 */
void Text_Token_Stream::pop_include() {
  Check(lines_.size() > 0);

  line_ = lines_.top();
  lines_.pop();
  buffer_ = buffers_.top();
  buffers_.pop();

  Require(check_class_invariants());
}

//------------------------------------------------------------------------------------------------//
Token Text_Token_Stream::scan_keyword() {
  Require(isalpha(peek_()) || peek_() == '_');

  string token_location = location_();
  /*char c =*/peek_();
  unsigned cc = 1;
  unsigned ci = 1;
  char c = peek_(ci);
  do {
    // Scan a C identifier.
    while (isalnum(c) || c == '_') {
      cc++;
      ci++;
      c = peek_(ci);
    }
    if (!no_nonbreaking_ws_) {
      // Replace any nonbreaking whitespace after the identifier with a single space, but ONLY if
      // the identifier is followed by another identifer.
      while (is_nb_whitespace(c)) {
        ci++;
        c = peek_(ci);
      }
      if (isalpha(c) || c == '_')
        cc++;
    }
  } while (isalpha(c) || c == '_');

  string text;
  text.reserve(cc);
  text += peek_(0);
  pop_char_();
  c = peek_();
  do {
    // Scan a C identifier.
    while (isalnum(c) || c == '_') {
      text += c;
      pop_char_();
      c = peek_();
    }
    if (!no_nonbreaking_ws_) {
      // Replace any nonbreaking whitespace after the identifier with a single space, but ONLY if
      // the identifier is followed by another identifer.
      while (is_nb_whitespace(c)) {
        pop_char_();
        c = peek_();
      }
      if (isalpha(c) || c == '_')
        text += ' ';
    }
  } while (isalpha(c) || c == '_');

  if (text == "end") {
    Ensure(check_class_invariants());
    return {END, give(token_location)};
  } else {
    Ensure(check_class_invariants());
    return {KEYWORD, give(text), give(token_location)};
  }
}

//------------------------------------------------------------------------------------------------//
Token Text_Token_Stream::scan_manifest_string() {
  Require(peek_() == '"');

  string token_location = location_();
  unsigned ci = 1;
  char c = peek_(ci);
  for (;;) {
    while (c != '"' && c != '\\' && c != '\n' && c != '\0') {
      ci++;
      c = peek_(ci);
    }
    if (c == '"')
      break;
    if (c == '\\') {
      ci += 2;
      c = peek_(ci);
    } else {
      if (c == '\0') {
        report_syntax_error(Token(EXIT, give(token_location)), "unexpected end of file; "
                                                               "did you forget a closing quote?");
      } else {
        Check(c == '\n');
        report_syntax_error(Token(EXIT, give(token_location)), "unexpected end of line; "
                                                               "did you forget a closing quote?");
      }
    }
  }
  ci++;

  string text;
  text.reserve(ci);
  text += peek_();
  pop_char_();
  c = peek_();
  for (;;) {
    while (c != '"' && c != '\\' && c != '\n' && c != '\0') {
      text += c;
      pop_char_();
      c = peek_();
    }
    if (c == '"')
      break;
    if (c == '\\') {
      text += c;
      pop_char_();
      c = pop_char_();
      text += c;
      c = peek_();
    } else {
      if (end_() || error_()) {
        report_syntax_error(Token(EXIT, give(token_location)), "unexpected end of file; "
                                                               "did you forget a closing quote?");
      } else {
        Check(c == '\n');
        report_syntax_error(Token(EXIT, give(token_location)), "unexpected end of line; "
                                                               "did you forget a closing quote?");
      }
    }
  }
  text += '"';
  pop_char_();

  //! \cond doxygen_ignore_block
  Ensure(check_class_invariants());
  //! \endcond

  return {STRING, give(text), give(token_location)};
}

} // end namespace rtt_parser

//------------------------------------------------------------------------------------------------//
// end of Text_Token_Stream.cc
//------------------------------------------------------------------------------------------------//
