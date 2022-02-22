//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   Token_Stream.cc
 * \author Kent G. Budge
 * \brief  Definitions of Token_Stream member functions.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

#include "Token_Stream.hh"
#include <cstring>

namespace rtt_parser {
using namespace std;

//------------------------------------------------------------------------------------------------//
Syntax_Error::Syntax_Error() : runtime_error("syntax error") {
  Ensure(!strcmp(what(), "syntax error"));
}

//------------------------------------------------------------------------------------------------//
/*!
 * This function returns the token at the cursor position and advance the
 * cursor. It will, if necessary, fill() the token buffer first.
 *
 * \return <code>old lookahead()</code>
 */
Token Token_Stream::shift() {
  if (deq.size() == 0) {
    deq.push_back(fill_());
  }

  Token Result;
  Result.swap(deq[0]);
  deq.pop_front();

  Ensure(check_class_invariants());
  // Ensure the cursor advances one place to the right, discarding the leftmost
  // token.
  return Result;
}

//------------------------------------------------------------------------------------------------//
/*!
 * This function looks ahead in the token stream without changing the cursor
 * position.  It will, if necessary, fill_() the token buffer first.  If the
 * requested position is at or past the end of the file, an EXIT token will be
 * returned.
 *
 * \param pos Number of tokens to look ahead, with 0 being the token at the
 * cursor position.
 *
 * \return The token at the specified position relative to the cursor.
 */
Token const &Token_Stream::lookahead(unsigned const pos) {
  while (deq.size() <= pos) {
    deq.push_back(fill_());
  }

  Ensure(check_class_invariants());
  return deq[pos];
}

//------------------------------------------------------------------------------------------------//
/*!
 * This function pushes the specified token onto the front of the token stream,
 * so that it is now the token in the lookahead(0) position.
 */
void Token_Stream::pushback(Token const &token) {
  deq.push_front(token);

  Ensure(check_class_invariants());
  Ensure(lookahead() == token);
}

//------------------------------------------------------------------------------------------------//
/*!
 * The default implementation of this function passes its message on to
 * Report_Error, then throws a Syntax_Error exception.
 *
 * A syntax error is a badly formed construct that requires explicit error
 * recovery (including resynchronization) by the parsing software.
 *
 * \param token
 * Token at which the error occurred.
 * \param message
 * The message to be passed to the user.
 *
 * \throw Syntax_Error This function never returns.  It always throws a
 * Syntax_Error exception to be handled by the parsing software.
 */
void Token_Stream::report_syntax_error(Token const &token, string const &message) {
  try {
    error_count_++;
    report(token, message);
  } catch (...) {
    // An error at this point really hoses us.  It means something went sour
    // with the reporting mechanism, and there probably isn't much we can do
    // about it.
    throw std::bad_exception();
  }

  Ensure(check_class_invariants());
  throw Syntax_Error();
}

//------------------------------------------------------------------------------------------------//
/*!
 * The default implementation of this function passes its message on to report,
 * then throws a Syntax_Error exception.
 *
 * A syntax error is a badly formed construct that requires explicit error
 * recovery (including resynchronization) by the parsing software.
 *
 * This versiona ssumes that the cursor is the location of the error.
 *
 * \param message
 * The message to be passed to the user.
 *
 * \throw Syntax_Error This function never returns.  It always throws a
 * Syntax_Error exception to be handled by the parsing software.
 */
void Token_Stream::report_syntax_error(string const &message) {
  try {
    error_count_++;
    report(message);
  } catch (...) {
    // An error at this point really hoses us.  It means something went sour
    // with the reporting mechanism, and there probably isn't much we can do
    // about it.
    throw std::bad_exception();
  }

  throw Syntax_Error();
}

//------------------------------------------------------------------------------------------------//
/*!
 * The default implementation of this function passes its message on to report,
 * then returns.
 *
 * A semantic error is a well-formed construct that has a bad value.  Because
 * the construct is well-formed, parsing may be able to continue after the error
 * is reported without any explicit recovery by the parsing software.
 *
 * \param token
 * Token at which the error occurred.
 * \param message
 * The message to be passed to the user.
 */
void Token_Stream::report_semantic_error(Token const &token, string const &message) {
  error_count_++;
  report(token, message);

  Ensure(check_class_invariants());
}

//------------------------------------------------------------------------------------------------//
/*!
 * The default implementation of this function passes its message on to report,
 * then returns.
 *
 * A semantic error is a well-formed construct that has a bad value.  Because
 * the construct is well-formed, parsing may be able to continue after the error
 * is reported without any explicit recovery by the parsing software.
 *
 * This version assumes that the cursor is the error location.
 *
 * \param message
 * The message to be passed to the user.
 */
void Token_Stream::report_semantic_error(string const &message) {
  error_count_++;
  report(message);

  Ensure(check_class_invariants());
}

//------------------------------------------------------------------------------------------------//
/*!
 *
 * The default implementation of this function passes its message on to report,
 * then returns.
 *
 * A semantic error is a well-formed construct that has a bad value.  Because
 * the construct is well-formed, parsing may be able to continue after the error
 * is reported without any explicit recovery by the parsing software.
 *
 * This version assumes that the cursor is the error location.
 *
 * \param message
 * The exception whose message is to be passed to the user.
 */
void Token_Stream::report_semantic_error(exception const &message) {
  error_count_++;
  report(message.what());

  Ensure(check_class_invariants());
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Reset the token stream.
 *
 * This function is normally called by its overriding version in children of
 * Token_Stream. It flushes the queues and resets the error count.
 */
void Token_Stream::rewind() {
  error_count_ = 0;
  deq.clear();

  Ensure(check_class_invariants());
}

} // namespace rtt_parser

//------------------------------------------------------------------------------------------------//
// end of Token_Stream.cc
//------------------------------------------------------------------------------------------------//
