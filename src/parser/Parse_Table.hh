//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   parser/Parse_Table.hh
 * \author Kent G. Budge
 * \brief  Definition of Keyword and Parse_Table.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_Parse_Table_HH
#define rtt_Parse_Table_HH

#include "Token_Stream.hh"
#include <cstring> // strcmp
#include <vector>

namespace rtt_parser {
//------------------------------------------------------------------------------------------------//
/*!
 * \brief Structure to describe a parser keyword.
 *
 * A Keyword describes a keyword in a Parse_Table.  It is a POD struct so that it can be initialized
 * using the low-level C++ initialization construct, e.g.,
 *
 * \code
 *   Keyword my_table[] = {{"FIRST",  Parse_First,  0, "TestModule"},
 *                         {"SECOND", Parse_Second, 0, "TestModule"}};
 * \endcode
 *
 * As a POD struct, Keyword can have no invariants.  However, Parse_Table imposes constraints on the
 * keywords it will accept for its keyword list.
 */
struct Keyword {
  /*! \brief The keyword moniker.
   *
   * The moniker should be a sequence of valid C++ identifiers separated by a single space.  For
   * example, \c "WORD", <tt>"First_Word Second_Word"</tt>, and <tt>"B1 T2 T3"</tt> are all valid
   * Keyword monikers.  Identifiers beginning with an underscore are permitted but may be reserved
   * for internal use by frameworks that uses the Parse_Table services. A Parse_Table attempts to
   * match input to the monikers in its Keyword table according to a set of rules stored in the
   * Parse_Table (q.v.)
   */
  char const *moniker = {nullptr};

  /*! \brief The keyword parsing function.
   *
   * When a Parse_Table finds a match to a moniker in its keyword table, the corresponding parse
   * function is called. The parse function may read additional tokens from the input Token_Stream,
   * such as a parameter value or an entire keyword block, before returning control to the
   * Parse_Table.
   *
   * \param stream
   * The token stream currently being parsed.
   *
   * \param index An integer argument that optionally allows a single parse function to handle a set
   * of related keywords. The
   */
  void (*func)(Token_Stream &stream, int index) = {nullptr};

  /*! \brief The index argument to the parse function.
   *
   * This is the index value that is passed to the parse function when the Parse_Table finds a match
   * to the keyword moniker. The parse function is not required to make any use of this argument,
   * but it can be convenient at times to use the same parse function for closely related keywords
   * and use the index argument to make the distinction.  For example, an enumerated or Boolean
   * option may be set using a single parse function that simply copies the index argument to the
   * option variable.
   */
  int index = {0};

  /*! Name of the module that supplied the keyword.
   *
   * This member is significant only if certain options are set in the Parse_Table.  It is used to
   * support parsing in frameworks that support self-registering modules.  Such modules often need
   * to add keywords to existing Parse_Tables, and the \c module member is useful for identifying
   * the module from which a particular keyword in a table originated.  This is particularly
   * important for diagnosing keyword clashes, where two modules have registered keywords with the
   * same moniker in the same Parse_Table.
   */
  char const *module = {nullptr};

  /*! Explanation of keyword.
   *
   * This optional member is a brief description of the keyword. If the parser sees a keyword in the
   * Token_Stream that it does not recognize, it will list the recognized keywords and, if a keyword
   * has a non-null description, it will print that description.
   */
  char const *description = {nullptr};

  Keyword() = default;

  Keyword(char const *moniker_in, void (*func_in)(Token_Stream &, int), int const index_in,
          char const *module_in, char const *description_in = nullptr)
      : moniker(moniker_in), func(func_in), index(index_in), module(module_in),
        description(description_in) {}
};

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Simple keyword-matching parse table
 *
 * A Parse_Table is a table of keywords and associated parsing functions.  It accepts tokens from a
 * Token_Stream, matching the keywords in the Token_Stream to its table and dispatching control to
 * the corresponding parsing functions. The parsing functions can take additional tokens, such as
 * parameter values, from the Token_Stream prior to returning control to the Parse_Table.
 *
 * How a Parse_Table determines whether an input keyword token matches a keyword moniker in its
 * Keyword table is determined by a set of flags.  By default, a Parse_Table is case-sensitive and
 * requires each identifier in the keyword token to exactly match the corresponding identifier in
 * the keyword moniker. However, a Parse_Table may be instructed to ignore case or to accept a
 * partial match on each identifier.  These options are described in more detail in connection with
 * the Set_Flags member function and the Keyword_Compare member functions.
 *
 * When the Parse_Table fails to match an input keyword to its keyword table, or when the input
 * token is not a keyword, END, EXIT, or ERROR, the parse table reports an error, then attempts
 * recovery by reading additional tokens until it encounters either a keyword it recognizes or an
 * END, EXIT, or ERROR. During this recovery, no additional errors are reported, to avoid swamping
 * the user with additional messages that are unlikely to be helpful.  If recovery is successful (a
 * keyword is recognized in the token stream) parsing resumes as normal.
 *
 * User parse routines may also encounter errors.  These may be reported to the Token_Stream through
 * either the Report_Syntax_Error or Report_Semantic_Error functions.  The former is used when there
 * is an error in the input syntax, as when a keyword is encountered when a numerical value is
 * expected.  In this case, an exception is thrown that is caught by the Parse_Table, which then
 * attempts error recovery as described above.  The Report_Semantic_Error indicates syntactically
 * correct input that is nonetheless unacceptable, as when a numerical value appears where it is
 * expected but violates some constraint, such as positivity.  The error is reported but input
 * processing is not interrupted.
 *
 * For parsers designed to receive input from the console, the default behavior just described is
 * not entirely appropriate. Console parsers are supported by treating any OTHER token whose text
 * character is a semicolon as an empty keyword.  By default, semicolons are treated as whitespace
 * and the parser will never see such a token.  A console stream can choose to convert an endline or
 * other terminator to a semicolon token to force processing.
 */
class Parse_Table {
public:
  // TYPEDEFS AND ENUMERATIONS

  //! Parser settings, which can be bitwise OR-ed together.
  enum {
    CASE_INSENSITIVE = 1U,         //!< Keyword match is case-insensitive
    PARTIAL_IDENTIFIER_MATCH = 2U, //!< Match incomplete identifiers
    ONCE = 4U                      //!< Parse single keyword, not a block
  };

  // CREATORS

  //! Create an empty Parse_Table.
  Parse_Table() : vec(0), flags_(0) {} // NOLINT

  //! Construct a parse table with the specified keywords.
  Parse_Table(Keyword const *table, size_t count, unsigned flags = 0);

  //! This class is meant to be heritable.
  virtual ~Parse_Table() = default;

  //! Disable copy construction
  Parse_Table(Parse_Table const &rhs) = delete;

  //! Disable move construction
  Parse_Table(Parse_Table &&rhs) noexcept = delete;

  //! Disable assignment
  Parse_Table &operator=(Parse_Table const &rhs) = delete;

  //! Disable move assignment
  Parse_Table &operator=(Parse_Table &&rhs) noexcept = delete;

  // MANIPULATORS

  //! Add keywords to the table.
  void add(Keyword const *table, size_t count) noexcept(false);

  //! Add the keywords from another Parse_Table
  void add(Parse_Table const &source) noexcept(false);

  //! Remove a keyword from the table.
  void remove(char const *moniker);

  //! Request a change in capacity.
  void reserve(std::vector<Keyword>::size_type n) { vec.reserve(n); }

  //! Set parser options.
  void set_flags(unsigned char f);

  // ACCESSORS

  //! Return the number of elements in the vector
  std::vector<Keyword>::size_type size() const { return vec.size(); }

  //! Return the current parser options.
  unsigned char get_flags() const { return flags_; }

  // SERVICES

  //! Parse a token stream.
  Token parse(Token_Stream &tokens) const;

  //! Parse a token stream.
  Token parseforkeyword(Token_Stream &tokens) const;

  //! Check the class invariants
  bool check_class_invariants() const;

private:
  // TYPEDEFS AND ENUMERATIONS

  //----------------------------------------------------------------------------------------------//
  /*!
   * \brief Ordering functor for Keyword
   *
   * Provides an ordering for Keyword compatible with STL sort and search routines.
   */
  class Keyword_Compare_ {
  public:
    explicit Keyword_Compare_(unsigned char flags);

    bool operator()(Keyword const &k1, Keyword const &k2) const;
    bool operator()(Keyword const &keyword, Token const &token) const noexcept;

    int kk_comparison(char const *m1, char const *m2) const;
    int kt_comparison(char const *m1, char const *m2) const noexcept;

  private:
    unsigned char flags_;
  };

  // IMPLEMENTATION

  //! Sort and check the table following the addition of new keywords
  void sort_table_() noexcept(false);

  // DATA

  std::vector<Keyword> vec;
  unsigned char flags_{0}; //!< Option flags for this parse table.
};

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Test equality of Keywords
 *
 * Tests equality of two Keyword structs. Note that this is strict equality, which differs from the
 * comparison offered by the Parse_Table::Keyword_Compare predicate class. It is used primarily for
 * checking assertions.
 *
 * \param a First keyword to be compared.
 * \param b Second keyword to be compared.
 *
 * \return <code>a.moniker == b.moniker &&a.func == b.func && a.index == b.index &&
 *	       a.module == b.module</code>
 */
inline bool operator==(Keyword const &a, Keyword const &b) {
  return strcmp(a.moniker, b.moniker) == 0 && a.func == b.func && a.index == b.index &&
         strcmp(a.module, b.module) == 0;
}

//------------------------------------------------------------------------------------------------//
//! Check whether a keyword is well-formed.
bool Is_Well_Formed_Keyword(Keyword const &key);

} // namespace rtt_parser

#endif // rtt_Parse_Table_HH

//------------------------------------------------------------------------------------------------//
// end of Parse_Table.hh
//------------------------------------------------------------------------------------------------//
