//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   Token_Stream.hh
 * \author Kent G. Budge
 * \brief  Definition of class Token_Stream.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_Token_Stream_HH
#define rtt_Token_Stream_HH

#include "Token.hh"
#include <deque>
#include <memory>

namespace rtt_parser {
//------------------------------------------------------------------------------------------------//
/*!
 * \brief Parser exception class
 *
 * This is an exception class for reporting syntax errors in simple parsers.
 */
class Syntax_Error : public std::runtime_error {
public:
  // CREATORS

  Syntax_Error();
};

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Abstract token stream for simple parsers
 *
 * Modern simulation codes require detailed problem specifications, which typically take the form of
 * a human-readable ASCII text input file.  The problem specification language used in such input
 * files has some similarities to a programming language, though it is typically simpler than a
 * high-level language like C++ or Java. The problem specification expressed in this problem
 * specification language must be scanned and parsed by the simulation code in order to load the
 * data structures and control parameters required to run the simulation.
 *
 * For example, if a fragment of a problem specification text takes the form
 *
 * <code>   boundary 1 gray source 1.0 keV  </code>
 *
 * then the simulation code must be able to recognize keywords like "boundary" or "gray source",
 * numerical values like "1" or "1.0", and expressions like "keV". It must also understand the
 * syntax of the problem specification language, which requires an integer value after a "boundary"
 * keyword, followed by a keyword identifying the kind of boundary ("gray source"), and so on.
 *
 * Modern input readers split the task of reading a problem specification text into scanning and
 * parsing. Scanning is the task of converting the raw text into a sequence of \b tokens, which
 * represent the keywords, numerical or string values, and other lowest-level constructs in the
 * problem specification language. This sequence or stream of tokens is then analyzed by a parser
 * that understands the syntax of the problem specification language and can extract the semantic
 * meaning of each part of the problem specification.
 *
 * A Token_Stream is an abstract representation of a stream of tokens that can be presented to a
 * Parse_Table or other parsing client. Each token is represented as a Token struct.  There is
 * unlimited lookahead and pushback capability, but no backtracking is implemented in this release.
 * In other words, the client may look as many tokens to the right of the cursor as he desires, and
 * he may push any number of tokens onto the stream at the cursor position; but when the cursor
 * advances (via the Shift method) the token it advances over is discarded forever.
 *
 * The actual scanner that does the conversion of raw text to tokens is provided by an
 * implementation of the \c fill function in a child class. This is usually the \c Text_Token_Stream
 * class, whose implementation of the \c fill function converts a stream of ASCII characters to a
 * stream of tokens, but one can imagine unconventional sources of tokens such as an HTML form on a
 * web interface.
 *
 * Because the token stream object is specific to a particular kind of user interface, we find it
 * convenient to allow the parser to report any syntax or semantic errors it discovers in the
 * problem specification to the token stream, which "knows" the best way to convey these to the
 * human client running the program. This style of error reporting also permits the token stream to
 * add additional information to the error message, indicating to the human client where the error
 * occurred. For example, a \c File_Token_Stream can tell the human client the line in the input
 * file where the error was detected.
 */

class Token_Stream {
public:
  // CREATORS

  virtual ~Token_Stream() = default;

  //! Copy constructor
  Token_Stream(Token_Stream const &rhs) = default;
  Token_Stream(Token_Stream &&rhs) = default;

  //! Copy assignment
  Token_Stream &operator=(Token_Stream const &rhs) = default;
  Token_Stream &operator=(Token_Stream &&rhs) = default;

  // MANIPULATORS

  //! Return the next token in the stream and advance the cursor.
  Token shift();

  /*!
   * \brief Look ahead at tokens.
   *
   * Lookahead references should remain valid until the referenced token is shifted.
   */
  Token const &lookahead(unsigned pos = 0);

  //! Insert a token into the stream at the cursor position.
  void pushback(Token const &token);

  //------------------------------------------------------------------------------------------------//
  /*!
   * \brief Reset the stream
   *
   * This function resets the token stream to some initial state defined by the
   * child class.
   */
  virtual void rewind() = 0;

  //! Report a syntax error to the user.
  [[noreturn]] void report_syntax_error(Token const &token, std::string const &message);

  //! Report a syntax error to the user.
  [[noreturn]] void report_syntax_error(std::string const &message);

  //! Report a semantic error to the user.
  virtual void report_semantic_error(Token const &token, std::string const &message);

  //! Report a semantic error to the user.
  virtual void report_semantic_error(std::string const &message);

  //! Report a semantic error to the user.
  virtual void report_semantic_error(std::exception const &message);

  //------------------------------------------------------------------------------------------------//
  /*!
   * \brief Report an error to the user.
   *
   * This function sends a message to the user in a stream-specific manner.
   *
   * \param token Token that triggered the error.
   *
   * \param message Message to be passed to the user.
   */
  virtual void report(Token const &token, std::string const &message) = 0;

  //------------------------------------------------------------------------------------------------//
  /*!
   * \brief Report an error to the user.
   *
   * This function sends a message to the user in a stream-specific manner. This variant assumes
   * that the cursor gives the correct location of the error.
   *
   * \param message
   * Message to be passed to the user.
   */
  virtual void report(std::string const &message) = 0;

  //------------------------------------------------------------------------------------------------//
  /*!
   * \brief Send a comment, without location information, to the user.
   *
   * This function sends a message to the user in a stream-specific manner. This differs from the
   * report() functions chiefly in that the message is not preceded by any location
   * information. This is useful for extended comments, e.g., listing available keywords when the
   * stream contains an unrecognized keyword.
   *
   * \param message
   * Message to be passed to the user.
   */
  virtual void comment(std::string const &message) = 0;

  //------------------------------------------------------------------------------------------------//
  /*! Check a syntax condition.
   *
   * By putting the check branch here, we improve coverage statistics for branch coverage. Making
   * the function inline keeps the cost of doing so negligible.
   *
   * \param condition Condition to be checked. If false, the message is delivered in a
   * stream-dependent manner, the stream's error counter is incremented, and a syntax exception is
   * thrown.
   *
   * \param message
   * Diagnostic message to be delivered if condition tests as false.
   */
  void check_syntax(bool const condition, char const *const message) {
    if (!condition)
      report_syntax_error(message);
  }

  //------------------------------------------------------------------------------------------------//
  /*! Check a semantic condition.
   *
   * By putting the check branch here, we improve coverage statistics for branch coverage. Making
   * the function inline keeps the cost of doing so negligible.
   *
   * \param condition Condition to be checked. If false, the message is delivered and the stream's
   * error counter is incremented.
   *
   * \param message Diagnostic message to be delivered if condition tests as false.
   */
  void check_semantics(bool const condition, char const *const message) {
    if (!condition)
      report_semantic_error(message);
  }

  // ACCESSORS

  //! Return the number of errors reported to the stream since it was last constructed or rewound.
  unsigned error_count() const noexcept { return error_count_; }

  //! Check that all class invariants are satisfied.
  bool check_class_invariants() const { return true; }

  // STATICS

protected:
  // IMPLEMENTATION

  //! Construct a Token_Stream.
  inline Token_Stream();

  //------------------------------------------------------------------------------------------------//
  /*!
   * \brief Add one or more tokens to the end of the lookahead buffer.
   *
   * This function is used by \c shift and \c lookahead to keep the token stream filled.
   *
   * Each call to the function scans the next token from the ultimate token source (such as a text
   * file or GUI) and returns it to the client. If no more tokens are available, the function must
   * return an EXIT token.
   *
   * \return Next token to put on the token buffer.
   */
  virtual Token fill_() = 0;

private:
  // DATA

  //! Number of errors reported to the stream since it was constructed or last rewound.
  unsigned error_count_ = {0};

private:
  // DATA
  std::deque<Token> deq;
};

//------------------------------------------------------------------------------------------------//
//! Construct a Token_Stream and place the cursor at the start of the stream.
inline Token_Stream::Token_Stream() : deq() {
  Ensure(check_class_invariants());
  Ensure(error_count() == 0);
}

} // namespace rtt_parser

#endif // rtt_Token_Stream_HH
//------------------------------------------------------------------------------------------------//
// end of Token_Stream.hh
//------------------------------------------------------------------------------------------------//
