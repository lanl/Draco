//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   Parse_Table.cc
 * \author Kent G. Budge
 * \date   Wed Jan 22 15:18:23 MST 2003
 * \brief  Definitions of member functions of class Parse_Table
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "Parse_Table.hh"
#include "ds++/Assert.hh"
#include <algorithm>
#include <cctype>
#include <sstream>

namespace rtt_parser {
using namespace std;

//------------------------------------------------------------------------------------------------//
/*!
 * \param table Pointer to an array of keywords.
 *
 * \param count Length of the array of keywords pointed to by \c table.
 * \param flags Initial value for parser flags.
 *
 * \throw invalid_argument If the keyword table is ill-formed or ambiguous.
 *
 * \note See documentation for \c Parse_Table::add for an explanation of the low-level argument
 *       list.
 */
Parse_Table::Parse_Table(Keyword const *const table, size_t const count, unsigned const flags)
    : vec(), flags_(static_cast<unsigned char>(flags)) {
  Require(count == 0 || table != nullptr);
  Require(count == 0 || std::find_if(table, table + count, Is_Well_Formed_Keyword));

  add(table, count);

  Ensure(check_class_invariants());
  Ensure(get_flags() == flags);
}

//------------------------------------------------------------------------------------------------//
/*!
 * \param table Array of keywords to be added to the table.
 * \param count Number of valid elements in the array of keywords.
 *
 * \throw invalid_argument If the keyword table is ill-formed or ambiguous.
 *
 * \note The argument list reflects the convenience of defining raw keyword tables as static C
 * arrays.  This justifies a low-level interface in place of, say, vector<Keyword>.
 */
void Parse_Table::add(Keyword const *const table, size_t const count) noexcept(false) {
  Require(count == 0 || table != nullptr);
  // Additional precondition checked in loop below

  // Preallocate storage.
  vec.reserve(vec.size() + count);

  // Add the new keywords.

  for (unsigned i = 0; i < count; i++) {
    Require(Is_Well_Formed_Keyword(table[i]));

    vec.push_back(table[i]);
  }

  sort_table_();

  Ensure(check_class_invariants());
}

//------------------------------------------------------------------------------------------------//
/*!
 * \param moniker Keyword to remove from the table.
 *
 * \throw invalid_argument If the keyword is not in the table.
 */
void Parse_Table::remove(char const *moniker) {
  // This is an order N operation as presently coded. N is never very large.
  for (auto i = vec.begin(); i != vec.end(); ++i) {
    if (!strcmp(i->moniker, moniker)) {
      vec.erase(i);
      Ensure(check_class_invariants());
      return;
    }
  }

  throw invalid_argument("keyword not found in Parse_Table::remove");
}

//------------------------------------------------------------------------------------------------//
/*!
 * \param source Parse_Table whose keywords are to be added to this Parse_Table.
 *
 * \throw invalid_argument If the keyword table is ill-formed or ambiguous.
 */
void Parse_Table::add(Parse_Table const &source) noexcept(false) {
  // Preallocate storage.
  vec.reserve(vec.size() + source.vec.size());

  // Add the new keywords.
  for (auto const &src : source.vec)
    vec.push_back(src);

  sort_table_();

  Ensure(check_class_invariants());
}

//------------------------------------------------------------------------------------------------//
/* private */
void Parse_Table::sort_table_() noexcept(false) // apparently std::sort can throw
{
  if (vec.size() == 0)
    return;

  // Sort the parse table, using a comparator predicate appropriate for the selected parser flags.

  Keyword_Compare_ const comp(flags_);
  std::sort(vec.begin(), vec.end(), comp);

  // Look for ambiguous keywords, and resolve the ambiguity, if possible.

  auto i = vec.begin();
  while (i + 1 != vec.end()) {
    Check(i->moniker != nullptr && (i + 1)->moniker != nullptr);
    if (!comp(i[0], i[1]))
    // kptr[i] and kptr[i+1] have the same moniker.
    {
      if (i->func == (i + 1)->func && i->index == (i + 1)->index) {
        // They have the same parse function and index.  No real ambiguity.  Delete the duplicate.

        // This can occur when there is diamond inheritance in a parse table hierarchy, e.g., this
        // parse table copies keywords from two other parse tables, which in turn copy keywords from
        // a fourth parse table.

        vec.erase(i + 1);
      } else {
        // The keywords are genuinely ambiguous. Throw an exception identifying the duplicate
        // keyword.
        using std::endl;
        using std::ostringstream;
        ostringstream err;
        err << "An ambiguous keyword was detected in a "
            << "Parse_Table:  " << i->moniker << endl
            << "Modules: ";
        if (i->module)
          err << i->module;
        else
          err << "<null>";
        err << ' ';
        if ((i + 1)->module)
          err << (i + 1)->module;
        else
          err << "<null>";
        err << endl;
        throw invalid_argument(err.str().c_str());
      }
    } else
    // kptr[i] and kptr[i+1] have different monikers. No ambiguity.
    {
      i++;
    }
  }
}
//------------------------------------------------------------------------------------------------//
/*!
 * Parse the stream of tokens until an END, EXIT, or ERROR token is reached.
 *
 * \param tokens The Token Stream from which to obtain the stream of tokens.
 *
 * \return The terminating token: either END, EXIT, or ERROR.
 *
 * \throw rtt_dsxx::assertion If the keyword table is ambiguous.
 */
Token Parse_Table::parse(Token_Stream &tokens) const {
  // The is_recovering flag is used during error recovery to suppress additional error messages.
  // This reduces the likelihood that a single error in a token stream will generate a large number
  // of error messages.

  bool is_recovering = false;

  // Create a comparator object that will be used to attempt to match keywords in the Token_Stream
  // to keywords in the Parse_Table.  This comparator object incorporates the current settings of
  // the Parse_Table, such as case sensitivity and partial matching options.

  Keyword_Compare_ const comp(flags_);

  // Now begin the process of pulling keywords off the input token stream, and attempting to match
  // these to the keyword table.

  for (;;) {
    Token token = tokens.shift();

    // The END, EXIT, and ERROR tokens are terminating tokens.  EXIT means the end of the token
    // stream has been reached.  END is used to flag the end of a nested parse, where the result of
    // matching a keyword in one parse table is to begin parsing keywords in a second parse table.
    // An END indicates that the second parse table should return control to the first parse table.
    // ERROR means that something went very wrong and we're probably hosed, but it allows some error
    // recovery from within a nested parse table.

    if (token.type() == END || token.type() == EXIT || token.type() == rtt_parser::ERROR) {
      return token;
    }

    // A Parse_Table assumes that every construct begins with a keyword.  This keyword is matched to
    // the keyword table, and if a match is found, control is directed to the associated parse
    // function, which can be written to accept just about any construct you wish.  However, by the
    // time return controls from a parse function, the token stream should be pointing either at a
    // terminating token or the next keyword.

    if (token.type() == KEYWORD) {
      // Attempt to match the keyword to the keyword table.  The following call returns an iterator
      // pointing to the first keyword in the table whose lexical ordering is greater than or equal
      // to the keyword token.  The lexical ordering is supplied by the comp object.

      vector<Keyword>::const_iterator const match =
          lower_bound(vec.begin(), vec.end(), token, comp);

      if (match == vec.end() || comp.kt_comparison(match->moniker, token.text().c_str()) != 0) {
        // The token was not lexically equal to anything in the keyword table.  In other words, the
        // keyword is unrecognized by the Parse_Table.  The error recovery procedure is to generate
        // a diagnostic, then pull additional tokens off the token stream (without generating
        // further diagnostics) until one is recognized as either a keyword or a terminating token.
        // We implement this behavior by setting the is_recovering flag when the first invalid token
        // is encountered, and resetting this flag as soon as a valid token is encountered.

        if (!is_recovering) {
          // We are not recovering from a previous error.  Generate a diagnostic, and flag that we
          // are now in error recovery mode.

          tokens.report_semantic_error(token, ": unrecognized keyword: " + token.text());

          // Give the user the possibilities.
          tokens.comment("Perhaps you meant one of:");
          for (auto const &i : vec) {
            tokens.comment(string("  ") + i.moniker);
            if (i.description != nullptr) {
              tokens.comment(string("    ") + i.description);
            }
          }

          is_recovering = true;
        }
        // else we are in recovery mode, and additional diagnostics are disabled until we see a
        // valid construct.
      } else {
        // We have a valid match.  However, depending on Parse_Table options, the match might be
        // ambiguous.  For example if the Parse_Table option to allow partial matches is active, the
        // keyword token may partially match more than one keyword in the keyword table.  Check for
        // an ambiguous match:

        if (match + 1 != vec.end() &&
            comp.kt_comparison(match[1].moniker, token.text().c_str()) == 0) {
          // The match is ambiguous.  This is diagnosed whether or not we are already in recovery
          // mode, but it does put us into recovery mode.

          tokens.report_semantic_error(token, "ambiguous keyword: " + token.text());
          is_recovering = true;
        } else {
          is_recovering = false;
          // We successfully processed something, so we are no longer in recovery mode.

          try {
            // Call the parse function associated with the keyword.
            match->func(tokens, match->index);

            if (flags_ & ONCE)
            // Quit after parsing a single keyword. This is useful for parse tables for selecting
            // one of a set of short options.
            {
              return {END, ""};
            }
          } catch (const Syntax_Error &) {
            // If the parse function detects a syntax error, and if it does not have its own error
            // recovery policy (or is unable to recover), it should call tokens.Report_Syntax_Error
            // which generates a diagnostic and throws a Syntax_Error exception. This puts the main
            // parser into recovery mode.

            is_recovering = true;
          }
        }
      }
    } else if (token.type() == OTHER && token.text() == ";") {
      // Treat a semicolon token as an empty keyword.  We are no longer in recovery mode, but we
      // don't actually do anything.

      is_recovering = false;
    } else {
      // The next token in the token stream is not a keyword, indicating a syntax error. Error
      // recovery consists of generating a diagnostic message, then continuing to pull tokens off
      // the token stream (without generating any further diagnostics) until one is recognized as
      // either a keyword or a terminating token.  We implement this behavior by setting the
      // is_recovering flag when the first invalid token is encountered, and resetting this flag as
      // soon as a valid token is encountered.

      if (!is_recovering) {
        // We are not recovering from a previous error.  Generate a diagnostic, and flag that we are
        // now in error recovery mode.

        std::ostringstream msg;
        msg << "expected a keyword, but saw " << token.text();
        tokens.report_semantic_error(token, msg.str());
        is_recovering = true;
      }
      // else we are in recovery mode, and additional diagnostics are disabled until we see a valid
      // construct.
    }
  }
}

//------------------------------------------------------------------------------------------------//
/*!
 * Parse the stream of tokens until a keyword is found or an END, EXIT, or ERROR token is reached.
 *
 * \param tokens The Token Stream from which to obtain the stream of tokens.
 *
 * \return The terminating token: either END, EXIT, or ERROR.
 *
 * \throw rtt_dsxx::assertion If the keyword table is ambiguous.
 */
Token Parse_Table::parseforkeyword(Token_Stream &tokens) const {
  // Create a comparator object that will be used to attempt to match keywords in the Token_Stream
  // to keywords in the Parse_Table.  This comparator object incorporates the current settings of
  // the Parse_Table, such as case sensitivity and partial matching options.

  Keyword_Compare_ const comp(flags_);

  // Now begin the process of pulling keywords off the input token stream, and attempting to match
  // these to the keyword table.

  for (;;) {
    Token token = tokens.shift();

    // The END, EXIT, and ERROR tokens are terminating tokens.  EXIT means the end of the token
    // stream has been reached.  END is used to flag the end of a nested parse, where the result of
    // matching a keyword in one parse table is to begin parsing keywords in a second parse table.
    // An END indicates that the second parse table should return control to the first parse table.
    // ERROR means that something went very wrong and we're probably hosed, but it allows some error
    // recovery from within a nested parse table.

    if (token.type() == END || token.type() == EXIT || token.type() == rtt_parser::ERROR) {
      return token;
    }

    // A Parse_Table assumes that every construct begins with a keyword.  This keyword is matched to
    // the keyword table, and if a match is found, control is directed to the associated parse
    // function, which can be written to accept just about any construct you wish.  However, by the
    // time return controls from a parse function, the token stream should be pointing either at a
    // terminating token or the next keyword.

    if (token.type() == KEYWORD) {
      // Attempt to match the keyword to the keyword table.  The following call returns an iterator
      // pointing to the first keyword in the table whose lexical ordering is greater than or equal
      // to the keyword token.  The lexical ordering is supplied by the comp object.

      vector<Keyword>::const_iterator const match =
          lower_bound(vec.begin(), vec.end(), token, comp);

      if (match == vec.end() || comp.kt_comparison(match->moniker, token.text().c_str()) != 0) {
        // The token was not lexically equal to anything in the keyword table.  In other words, the
        // keyword is unrecognized by the Parse_Table.  The error recovery procedure is to generate
        // a diagnostic, then pull additional tokens off the token stream (without generating
        // further diagnostics) until one is recognized as either a keyword or a terminating token.
        // We implement this behavior by setting the is_recovering flag when the first invalid token
        // is encountered, and resetting this flag as soon as a valid token is encountered.
      } else {
        // We have a valid match.  However, depending on Parse_Table options, the match might be
        // ambiguous.  For example if the Parse_Table option to allow partial matches is active, the
        // keyword token may partially match more than one keyword in the keyword table.  Check for
        // an ambiguous match:

        if (match + 1 != vec.end() &&
            comp.kt_comparison(match[1].moniker, token.text().c_str()) == 0) {
          // The match is ambiguous.  This is diagnosed whether or not we are already in recovery
          // mode, but it does put us into recovery mode.

          tokens.report_semantic_error(token, "ambiguous keyword: " + token.text());
        } else {
          // We successfully processed something, so we are no longer in recovery mode.

          try {
            // Call the parse function associated with the keyword.
            match->func(tokens, match->index);

            if (flags_ & ONCE)
            // Quit after parsing a single keyword. This is useful for parse tables for selecting
            // one of a set of short options.
            {
              return {END, ""};
            }
          } catch (const Syntax_Error &) {
            // If the parse function detects a syntax error, and if it does not have its own error
            // recovery policy (or is unable to recover), it should call tokens.Report_Syntax_Error
            // which generates a diagnostic and throws a Syntax_Error exception. This puts the main
            // parser into recovery mode.
            tokens.report_semantic_error(token, "syntax error: " + token.text());
          }
        }
      }
    }
  }
}

//------------------------------------------------------------------------------------------------//
/*!
 * \param f Flags to be set, ORed together.
 */
void Parse_Table::set_flags(unsigned char const f) {
  flags_ = f;

  add(nullptr, 0U);
  // The keyword list needs to be sorted and checked.  For example, if the options are changed so
  // that a previously case-sensitive Parse_Table is no longer case-sensitive, then the ordering
  // changes, and previously unambiguous keywords may become ambiguous.

  Ensure(check_class_invariants());
  Ensure(get_flags() == f);
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Constructor for comparison predicate for sorting keyword tables.
 *
 * The predicate is used by a Parse_Table to sort its keyword list using std::sort.
 *
 * \param flags The flags controlling this comparator's operations.
 */
Parse_Table::Keyword_Compare_::Keyword_Compare_(unsigned char const flags) : flags_(flags) {}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Comparison function for sorting keyword tables.
 *
 * This function is used by a Parse_Table to sort its keyword list using std::sort.
 *
 * If no option flags are set, monikers will test equal if they are identical.
 *
 * If the CASE_INSENSITIVE option is set, monikers will test equal if they are identical when
 * converted to uppercase.
 *
 * Note that PARTIAL_IDENTIFIER_MATCH is ignored when comparing monikers from two keywords. It is
 * relevant only when comparing monikers with input tokens.
 *
 * A valid Parse_Table may not contain any keywords that test equal.
 *
 * \param k1 The first Keyword to be compared.
 *
 * \param k2 The second Keyword to be compared.
 *
 * \return <CODE>kk_comparison(k1.moniker, k2.moniker)<0 </CODE>
 */
bool Parse_Table::Keyword_Compare_::operator()(Keyword const &k1, Keyword const &k2) const {
  Require(k1.moniker != nullptr);
  Require(k2.moniker != nullptr);

  return kk_comparison(k1.moniker, k2.moniker) < 0;
}

int Parse_Table::Keyword_Compare_::kk_comparison(char const *m1, char const *m2) const {
  using namespace std;

  Require(m1 != nullptr);
  Require(m2 != nullptr);

  if (flags_ & CASE_INSENSITIVE) {
    while (*m1 != '\0' && *m2 != '\0') {
      char c1 = *m1++;
      char c2 = *m2++;
      if (islower(c1))
        c1 = static_cast<char>(::toupper(c1));
      if (islower(c2))
        c2 = static_cast<char>(::toupper(c2));
      if (c1 < c2)
        return -1;
      if (c1 > c2)
        return 1;
    }
    char c1 = *m1;
    char c2 = *m2;
    if (c1 < c2)
      return -1;
    if (c1 > c2)
      return 1;
    return 0;
  } else {
    return strcmp(m1, m2);
  }
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Comparison function for finding token match in keyword table.
 *
 * This function is used by a Parse_Table to match keywords to identifier tokens using
 * std::lower_bound.
 *
 * If no option flags are set, the match must be exact.
 *
 * If CASE_INSENSITIVE is set, the match must be exact after conversion to uppercase.
 *
 * If PARTIAL_IDENTIFIER_MATCH is set, each identifier in the token must be a prefix of the
 * corresponding identifier in the keyword, after conversion to uppercase if CASE_INSENSITIVE is
 * also set.  For example, the token "ABD" matches the keyword "ABDEF" but not the keyword "AB".
 *
 * \param k1 The Keyword to be compared.
 * \param k2 The token to be compared.
 *
 * \return <CODE>comparison(keyword.moniker, token.text().c_str())<0 </CODE>
 */
bool Parse_Table::Keyword_Compare_::operator()(Keyword const &k1, Token const &k2) const noexcept {
  Require(k1.moniker != nullptr);

  return kt_comparison(k1.moniker, k2.text().c_str()) < 0;
}

int Parse_Table::Keyword_Compare_::kt_comparison(char const *m1, char const *m2) const noexcept {
  using namespace std;

  Require(m1 != nullptr);
  Require(m2 != nullptr);

  if (flags_ & PARTIAL_IDENTIFIER_MATCH) {
    while (*m1 != '\0' && *m2 != '\0') {
      while (isalnum(*m1) && isalnum(*m2)) {
        char c1 = *m1++;
        char c2 = *m2++;
        if (flags_ & CASE_INSENSITIVE) {
          if (islower(c1))
            c1 = static_cast<char>(::toupper(c1));
          if (islower(c2))
            c2 = static_cast<char>(::toupper(c2));
        }
        if (c1 < c2)
          return -1;
        if (c1 > c2)
          return 1;
      }
      if (*m1 == ' ' && *m2 != ' ')
        return 1;
      while (isalnum(*m1))
        m1++;
      while (isalnum(*m2))
        m2++;
      if (*m1 == ' ')
        ++m1;
      if (*m2 == ' ')
        ++m2;
    }
    if (*m1 == '\0' && *m2 != '\0')
      return -1;
    if (*m1 != '\0' && *m2 == '\0')
      return 1;
    Check(*m1 == '\0' && *m2 == '\0');
    return 0;
  } else {
    if (flags_ & CASE_INSENSITIVE) {
      while (*m1 != '\0' && *m2 != '\0') {
        char c1 = *m1++;
        char c2 = *m2++;
        if (islower(c1))
          c1 = static_cast<char>(::toupper(c1));
        if (islower(c2))
          c2 = static_cast<char>(::toupper(c2));
        if (c1 < c2)
          return -1;
        if (c1 > c2)
          return 1;
      }
      char c1 = *m1++;
      char c2 = *m2++;
      if (islower(c1))
        c1 = static_cast<char>(::toupper(c1));
      if (islower(c2))
        c2 = static_cast<char>(::toupper(c2));
      if (c1 < c2)
        return -1;
      if (c1 > c2)
        return 1;
      return 0;
    } else {
      return strcmp(m1, m2);
    }
  }
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Check whether a keyword satisfies the requirements for use in a Parse_Table.
 *
 * \param key Keyword to be checked.
 *
 * \return \c false unless all the following conditions are met: <ul><li>\c key.moniker must point
 * to a null-terminated string consisting of one or more valid C++ identifiers separated by single
 * spaces.</li> <li>\c key.func must point to a parsing function.</li></ul>
 */
bool Is_Well_Formed_Keyword(Keyword const &key) {
  using namespace std;

  if (key.moniker == nullptr || key.func == nullptr)
    return false;
  char const *cptr = key.moniker;
  for (;;) {
    // Must be at the start of a C identifier, which begins with an alphabetic character or an
    // underscore.
    if (*cptr != '_' && !isalpha(*cptr))
      return false;

    // The subsequent characters in a C identifier must be alphanumeric characters or underscores.
    while (*cptr == '_' || isalnum(*cptr))
      cptr++;

    // If the identifier is followed by a null, we're finished scanning a valid keyword.
    if (*cptr == '\0')
      return true;

    // Otherwise, if the next character is not a space, it's not a valid keyword.
    if (*cptr != ' ')
      return false;

    // Skip over the space. cptr should now point to the start of the next C identifier, if this is
    // a valid keyword.
    cptr++;
  }
}

//------------------------------------------------------------------------------------------------//
bool Parse_Table::check_class_invariants() const {
  // The keyword table must be well-formed, sorted, and unambiguous.

  Keyword_Compare_ const comparator(flags_);
  for (auto i = vec.begin(); i != vec.end(); ++i) {
    if (!Is_Well_Formed_Keyword(i[0]))
      return false;
    if (i + 1 != vec.end()) {
      if (comparator.kk_comparison(i[0].moniker, i[1].moniker) >= 0)
        return false;
    }
  }
  return true;
}

} // namespace rtt_parser

//------------------------------------------------------------------------------------------------//
// end of Parse_Table.cc
//------------------------------------------------------------------------------------------------//
