//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   Class_Parser.i.hh
 * \brief  Definitions of member functions of template Class_Parser
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_Class_Parser_i_hh
#define rtt_Class_Parser_i_hh

#include <sstream>

namespace rtt_parser {

//------------------------------------------------------------------------------------------------//
/*!
 * \param child Reference to the complete child object for which this base is being constructed.
 *
 * \param raw_table Pointer to an array of keywords.
 *
 * \param[in] count Length of the array of keywords pointed to by \c table.
 *
 * \throw invalid_argument If the keyword table is ill-formed or ambiguous.
 *
 * \note See documentation for \c Parse_Table::add for an explanation of the low-level argument
 *       list.
 */
template <typename Class, bool once, bool allow_exit>
Class_Parser_Base<Class, once, allow_exit>::Class_Parser_Base(
    Child &child, Class_Parser_Keyword<Class> const *raw_table, unsigned const count)
    : child_(child) {
  Require(count == 0 || std::find_if(raw_table, raw_table + count, Is_Well_Formed_Keyword<Class>));

  static bool first_time = true;
  if (first_time) {
    add(raw_table, count);
    first_time = false;
  }

  Ensure(check_class_invariants());
}

//------------------------------------------------------------------------------------------------//
/*!
 * \param[in] table Array of keywords to be added to the table.
 * \param[in] count Number of valid elements in the array of keywords.
 *
 * \throw invalid_argument If the keyword table is ill-formed or ambiguous.
 *
 * \note The argument list reflects the convenience of defining raw keyword tables as static C
 *       arrays.  This justifies a low-level interface in place of, say, vector<Keyword>.
 */
template <typename Class, bool once, bool allow_exit>
void Class_Parser_Base<Class, once, allow_exit>::add(Keyword const *const table,
                                                     size_t const count) noexcept(false) {
  Require(count == 0 || table != nullptr);
  // Additional precondition checked in loop below

  if (count > 0) {
    // Preallocate storage.
    reserve(size() + count);

    // Add the new keywords.

    for (unsigned i = 0; i < count; i++) {
      Require(Is_Well_Formed_Keyword(table[i]));

      table_.push_back(table[i]);
    }

    sort_table_();
  }

  Ensure(check_class_invariants());
}

//------------------------------------------------------------------------------------------------//
/*!
 * This function is provided to support inheritance of class parsers.
 *
 * \param table Vector of keywords to be added to the table.
 *
 * \throw invalid_argument If the keyword table is ill-formed or ambiguous.
 */
template <typename Class, bool once, bool allow_exit>
template <typename Base>
void Class_Parser_Base<Class, once, allow_exit>::add(
    std::vector<Class_Parser_Keyword<Base>> const &table) noexcept(false) {
  unsigned const count = table.size();
  // Additional precondition checked in loop below

  if (count > 0) {
    // Preallocate storage.
    reserve(size() + count);

    // Add the new keywords.

    for (unsigned i = 0; i < count; i++) {
      auto const &b = table[i];
      Keyword d(b.moniker, b.func, b.index, b.module, b.description);
      Require(Is_Well_Formed_Keyword(d));
      table_.push_back(d);
    }

    sort_table_();
  }

  Ensure(check_class_invariants());
}

//------------------------------------------------------------------------------------------------//
/* private */
template <typename Class, bool once, bool allow_exit>
void Class_Parser_Base<Class, once, allow_exit>::sort_table_() noexcept(
    false) // apparently std::sort can throw
{
  Require(table_.size() > 0);

  // Sort the parse table, using a comparator predicate appropriate for the selected parser flags.

  Keyword_Compare_ const comp;
  std::sort(table_.begin(), table_.end(), comp);

  // Look for ambiguous keywords, and resolve the ambiguity, if possible.

  auto i = table_.begin();
  while (i + 1 != table_.end()) {
    Check(i->moniker != nullptr && (i + 1)->moniker != nullptr);
    if (!comp(i[0], i[1]))
    // kptr[i] and kptr[i+1] have the same moniker.
    {
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
      throw std::invalid_argument(err.str().c_str());
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
template <typename Class, bool once, bool allow_exit>
std::shared_ptr<Class> Class_Parser_Base<Class, once, allow_exit>::parse(Token_Stream &tokens) {

  using std::vector;

  // Save the old error count, so we can distinguish fresh errors within this class keyword block
  // from previous errors.
  unsigned const old_error_count = tokens.error_count();

  // The is_recovering flag is used during error recovery to suppress additional error messages.
  // This reduces the likelihood that a single error in a token stream will generate a large number
  // of error messages.

  bool is_recovering = false;

  // Create a comparator object that will be used to attempt to match keywords in the Token_Stream
  // to keywords in the Parse_Table.

  Keyword_Compare_ const comp;

  // Now begin the process of pulling keywords off the input token stream, and attempting to match
  // these to the keyword table.

  for (;;) {
    Token const token = tokens.shift();

    // The END, EXIT, and ERROR tokens are terminating tokens.  EXIT means the end of the token
    // stream has been reached.  END is used to flag the end of a nested parse, where the result of
    // matching a keyword in one parse table is to begin parsing keywords in a second parse table.
    // An END indicates that the second parse table should return control to the first parse table.
    // ERROR means that something went very wrong and we're probably hosed, but it allows some error
    // recovery from within a nested parse table.

    if (token.type() == END || token.type() == EXIT || token.type() == rtt_parser::ERROR) {
      if (token.type() == END || (allow_exit && token.type() == EXIT))
      // A class keyword block is expected to end with an END or (if allow_exit is true) an EXIT.
      {
        check_completeness(tokens);

        if (tokens.error_count() == old_error_count) {
          // No fresh errors in the class keyword block.  Create the object.
          return std::shared_ptr<Class>(create_object());
        } else {
          // there were errors in the keyword block. Don't try to create a class object.  Return the
          // null pointer.
          return nullptr;
        }
      } else {
        tokens.report_syntax_error("missing 'end'?");
        // never returns
      }
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

      auto const match = lower_bound(table_.begin(), table_.end(), token, comp);

      if (match == table_.end() || strcmp(match->moniker, token.text().c_str()) != 0) {
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
          for (auto const &i : table_) {
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
        // We have a valid match.

        is_recovering = false;
        // We successfully processed something, so we are no longer in recovery mode.

        try {
          // Call the parse function associated with the keyword.
          (child_.*(match->func))(tokens, match->index);

          if (once)
          // Quit after parsing a single keyword. This is useful for parse tables for selecting
          // one of a set of short options.
          {
            check_completeness(tokens);

            if (tokens.error_count() == old_error_count) {
              // No fresh errors in the class keyword block.  Create the object.
              return std::shared_ptr<Class>(create_object());
            } else {
              // there were errors in the keyword block. Don't try to create a class object.  Return
              // the null pointer.
              return nullptr;
            }
          }
        } catch (const Syntax_Error &) {
          // If the parse function detects a syntax error, and if it does not have its own error
          // recovery policy (or is unable to recover), it should call tokens.Report_Syntax_Error
          // which generates a diagnostic and throws a Syntax_Error exception. This puts the main
          // parser into recovery mode.

          is_recovering = true;
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
 * \brief Comparison function for sorting keyword tables.
 *
 * This function is used by a Parse_Table to sort its keyword list using std::sort.
 *
 * A valid Parse_Table may not contain any keywords that test equal.
 *
 * \param k1 The first Keyword to be compared.
 *
 * \param k2 The second Keyword to be compared.
 *
 * \return <CODE>strcmp(k1.moniker, k2.moniker)<0 </CODE>
 */
template <typename Class, bool once, bool allow_exit>
bool Class_Parser_Base<Class, once, allow_exit>::Keyword_Compare_::operator()(
    Keyword const &k1, Keyword const &k2) const {
  Require(k1.moniker != nullptr);
  Require(k2.moniker != nullptr);
  return strcmp(k1.moniker, k2.moniker) < 0;
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Comparison function for finding token match in keyword table.
 *
 * This function is used by a Parse_Table to match keywords to identifier tokens using
 * std::lower_bound.
 *
 * \param k1 The Keyword to be compared.
 * \param k2 The token to be compared.
 *
 * \return <CODE>strcmp(keyword.moniker, token.text().c_str())<0 </CODE>
 */
template <typename Class, bool once, bool allow_exit>
bool Class_Parser_Base<Class, once, allow_exit>::Keyword_Compare_::operator()(
    Keyword const &k1, Token const &k2) const noexcept {
  Require(k1.moniker != nullptr);
  return strcmp(k1.moniker, k2.text().c_str()) < 0;
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Check whether a keyword satisfies the requirements for use in a Parse_Table.
 *
 * \param key Keyword to be checked.
 *
 * \return \c false unless all the following conditions are met: (1) \c key.moniker must point to a
 *         null-terminated string consisting of one  or more valid C++ identifiers separated by
 *         single spaces. (2) \c key.func must point to a parsing function.
 */
template <typename Class> bool Is_Well_Formed_Keyword(Class_Parser_Keyword<Class> const &key) {
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
template <typename Class, bool once, bool allow_exit>
bool Class_Parser_Base<Class, once, allow_exit>::check_class_invariants() const {
  // The keyword table must be well-formed, sorted, and unambiguous.

  for (auto i = table_.begin(); i != table_.end(); ++i) {
    if (!Is_Well_Formed_Keyword(i[0]))
      return false;
    if (i + 1 != table_.end()) {
      if (strcmp(i[0].moniker, i[1].moniker) >= 0)
        return false;
    }
  }
  return true;
}

//------------------------------------------------------------------------------------------------//
template <typename Class, bool once, bool allow_exit>
std::vector<typename Class_Parser_Base<Class, once, allow_exit>::Keyword>
    Class_Parser_Base<Class, once, allow_exit>::table_;

} // namespace rtt_parser

#endif // rtt_Class_Parser_i_hh

//------------------------------------------------------------------------------------------------//
// end of Class_Parser.i.hh
//------------------------------------------------------------------------------------------------//
