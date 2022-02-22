//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   parser/Class_Parser.hh
 * \brief  Definition of Class_Parser_Keyword and Class_Parser.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_Class_Parser_HH
#define rtt_Class_Parser_HH

#include "Token_Stream.hh"
#include <cstring> // strcmp
#include <memory>
#include <vector>

namespace rtt_parser {

//------------------------------------------------------------------------------------------------//
/*! Template for parser table classes
 *
 * A Class_Parser is conceptually a table of keywords and associated parsing functions. When invoked
 * (via a call to the Class_Parser::parse() member function), it pulls tokens from a Token_Stream,
 * and attempts to match the tokens against a table of keywords.  Each keyword has a parsing
 * function and integer flag associate with it, and when a keyword matches an input token, the
 * parsing function is called with the original Token_Stream and the integer flag as arguments. The
 * parsing functions can pull additional tokens from the Token_Stream, such as parameter values,
 * prior to returning control to the Class_Parser. The machinery to do all this is normally provided
 * by deriving the Class_Parser from Class_Parser_Base.
 *
 * Each parse function is a member function of the Class_Parser, and the values it parses are
 * normally stored in the Class_Parser for later use to construct an object of the class for which
 * the Class_Parser is specialized. A pointer to this object is returned by parse().
 *
 * When the Class_Parse_Table_Base fails to match an input keyword to its keyword table, or when the
 * input token is not a keyword, END, EXIT, or ERROR, the parse table reports an error, then
 * attempts recovery by reading additional tokens until it encounters either a keyword it recognizes
 * or an END, EXIT, or ERROR.  During this recovery, no additional errors are reported, to avoid
 * swamping the user with additional messages that are unlikely to be helpful.  If recovery is
 * successful (a keyword is recognized in the token stream) parsing resumes to diagnose any
 * additional errors.
 *
 * Keyword parse routines may also encounter errors.  These may be reported to the Token_Stream
 * through either the report_syntax_error or report_semantic_error functions.  The former is used
 * when there is an error in the input syntax, as when a keyword is encountered when a numerical
 * value is expected.  In this case, an exception is thrown that is caught by the Class_Parser_Base,
 * which then attempts error recovery as described above.  The report_semantic_error indicates
 * syntactically correct input that is nonetheless unacceptable, as when a numerical value appears
 * where it is expected but violates some constraint, such as positivity.  The error is reported but
 * input processing is not interrupted. In effect, this is an error with immediate recovery.
 *
 * No general implementation is provided. The developer wishing to make his class parsable must
 * specialize this template for his class. This should almost always be derived from the
 * Class_Parser_Base template, which provides the necessary boilerplate. The final specialization
 * must define the following members:
 *
 *   Class_Parser()
 *
 * Default constructor. The constructor should set all data members to sentinel values indicating
 * that the corresponding parameters have not been parsed yet.  Where possible, the sentinel value
 * should be one that doesn't make sense for the parameter in question, such as 0 for a dimension
 * that must be nonzero, or a negative value for a parameter that must be positive.  This ensures
 * that the sentinel value can be distinguished from any value a user might specify in an input
 * "deck." If the parameter can plausibly have any value, then a bool sentinel flag should accompany
 * the parameter and be initialized in the constructor to false (for not yet parsed).
 *
 * The constructor should also carry out any other preparations that are necessary before parsing.
 * The constructor may optionally have a context argument to provide any contextual information
 * required to know how to parse the "deck".
 *
 *   void check_completeness(Token_Stream &);
 *
 * This function is called to check that all required specifications were found in the parsed token
 * stream.  This function must call \c Token_Stream::report_semantic_error at least once if a
 * required specification is missing.
 *
 * If default values are permitted for some parameters (which is not recommended), then they should
 * be applied here.
 *
 *   Class *create_object();
 *
 * Create the object from the parsed fields.  This function should have no preconditions that are
 * not guaranteed by a preceding successful call to check_completeness. create_object will be called
 * only if no errors were detected in the input token stream.
 *
 * In addition to these mandatory members, any nontrivial Class_Parser will also define data members
 * to hold the parsed specifications, and parse functions of the form
 *
 * void parse_parameter(Token_Stream &, int);
 *
 * A nontrivial Class_Parser should define a member
 *
 * \code
 * static Class_Parser_Keyword<Class> const raw_table[];
 * \endcode
 *
 * to hold a raw table of keywords that is passed to the constructor for Class_Parser_Base in the
 * constructor for Class_Parser. The Class_Parser_Base will inspect this table for errors and make a
 * sorted copy for fast token matching.
 *
 * The specialization will thus look something like
 *
 * Class__parser.hh:
 *
 * \code
 * template<>
 * class Class_Parser<Class> : public Class_Parser_Base<Class> {
 *    public:
 *      explicit Class_Parser(Context_Type debug_context); // context argument is optional
 *
 *      void check_completeness(Token_Stream &tokens);
 *
 *      Class *create_object();
 *
 *    protected:
 *      // Data to be parsed which will be used to construct the class, for example:
 *
 *      bool flag;
 *      double constant;
 *
 *      // Optional context information
 *      Context_Type debug_context;
 *
 *    private:
 *      // Parse functions
 *      friend Class_Parser_Base<Class>;
 *      void parse_flag(Token_Stream &, int);
 *      void parse_constant(Token_Stream &, int);
 *
 *      // The raw keyword table.
 *      Class_Parser_Keyword<Class> const raw_table[];
 * };
 * \endcode
 *
 *
 * Class__parser.cc:
 *
 * \code
 * include "Class__parser.hh"
 *
 * void Class_Parser<Class>::parse_flag(Token_Stream &tokens, int) {
 *    flag = parse_bool(tokens);
 * }
 *
 * void Class_Parser<Class>::parse_constant(Token_Stream &tokens, int) {
 *    constant = parse_real(tokens);
 * }
 *
 * template<>
 * Class_Parser_Keyword<Class> const Class_Parser_Base<Class>::raw_table[] = {
 *     {"block", &Class_Parser<Class>::parse_block, 0, ""},
 *     {"solver", &Class_Parser<Class>::parse_solver, 0, ""},
 *   };
 * // Note that this must be visible before the constructor so that the constructor can compute its
 * // length correctly.
 *
 * Class_Parser<Class>::Class_Parser(Context_Type debug_word)
 * : Class_Parser_Base<Class>(raw_table, sizeof(raw_table)/sizeof(raw_table[0]),
 *   debug_context(debug_word), flag(false), constant(0.0) {
 * }
 *
 * void Class_Parse_Table<Class>::check_completeness(Token_Stream &tokens) {
 *   ... check the parsed data ...
 * }
 *
 * Class *Class_Parse_Table<Class>::create_object() {
 *   return make_shared<Class>(debug_context, flag, constant);
 * }
 *
 * template <>
 * shared_ptr<Class> parse_class<Class>(Token_Stream &tokens, Debug_Context debug_context) {
 *    Class_Parser<Class> parser(debug_context);
 *
 *    return parser.parse(tokens);
 * }
 * \endcode
 *
 * This introduces all the "boilerplate" and lets the developer focus on the data required for the
 * constructor for Class, the parse functions needed to parse this data, and the check and construct
 * functions.
 *
 * Inheritance is supported, except for diamond inheritance (virtual base classes).
 */

template <typename Class> class Class_Parser;

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Structure to describe a parser keyword.
 *
 * A Class_Parser_Keyword describes a keyword in a Class_Parser.  It is a POD struct so that it can
 * be initialized using the low-level C++ initialization construct, e.g.,
 *
 * \code
 *   Class_Parser_Keyword<My_Class> my_table[] = {
 *      {"FIRST",  Parse_First,  0, "TestModule"},
 *      {"SECOND", Parse_Second, 0, "TestModule"}
 *   };
 * \endcode
 *
 * As a POD struct, Class_Parser_Keyword can have no invariants.  However, Class_Parser imposes
 * constraints on the keywords it will accept for its keyword list.
 */
template <typename Class> struct Class_Parser_Keyword {
  /*! \brief The keyword moniker.
   *
   * The moniker should be a sequence of valid C++ identifiers separated by a single space.  For
   * example, \c "WORD", \c "First_Word Second_Word", and \code "B1 T2 T3" \endcode are all valid
   * Keyword monikers.  Identifiers beginning with an underscore are permitted but may be reserved
   * for internal use by frameworks that uses the Class_Parser services. A Class_Parser attempts to
   * match input to the monikers in its Class_Parser_Keyword table according to a set of rules
   * stored in the Class_Parser (q.v.)
   */
  char const *moniker;

  /*! \brief The keyword parsing function.
   *
   * When a Class_Parser finds a match to a moniker in its keyword table, the corresponding parse
   * function is called. The parse function may read additional tokens from the input
   * Token_Stream, such as a parameter value or an entire keyword block, before returning control
   * to the Class_Parser.
   *
   * \param stream
   * The token stream currently being parsed.
   *
   * \param index An integer argument that optionally allows a single parse function to handle a
   * set of related keywords. The
   */
  void (Class_Parser<Class>::*func)(Token_Stream &stream, int index);

  /*! \brief The index argument to the parse function.
   *
   * This is the index value that is passed to the parse function when the
   * Parse_Table finds a match to the keyword moniker. The parse function is
   * not required to make any use of this argument, but it can be convenient
   * at times to use the same parse function for closely related keywords
   * and use the index argument to make the distinction.  For example, an
   * enumerated or Boolean option may be set using a single parse function
   * that simply copies the index argument to the option variable.
   */
  int index;

  /*! Name of the module that supplied the keyword.
   *
   * This member supports diagnostics for inherited Class_Parsers.
   * It identifies which class in the inheritance hierarchy is associated
   * with a particular keyword.
   */
  char const *module;

  /*! Explanation of keyword.
   *
   * This optional member is a brief description of the keyword. If the
   * parser sees a keyword in the Token_Stream that it does not recognize,
   * it will list likely matches ("Did you mean ..."). If a keyword has a
   * non-null description, it will print that description.
   */
  char const *description;

  Class_Parser_Keyword(char const *moniker, void (Class_Parser<Class>::*func)(Token_Stream &, int),
                       int const index, char const *module, char const *description = nullptr)
      : moniker(moniker), func(func), index(index), module(module), description(description) {}
};

//-------------------------------------------------------------------------//
/*!
 * \brief Boilerplate base for Class_Parser specializations
 *
 * Class_Parse_Base provides the boilerplate underlying almost every
 * Class_Parser specialization.
 */

template <typename Class, bool once = false, bool allow_exit = false> class Class_Parser_Base {
public:
  // TYPEDEFS AND ENUMERATIONS

  typedef Class_Parser<Class> Child;
  typedef Class_Parser_Keyword<Class> Keyword;
  typedef void (Child::*Parse_Function)(Token_Stream &, int);

  // CREATORS

  //! Construct a parse table with no keywords. Useful as a placeholder.
  Class_Parser_Base(Child &child) : child_(child) {}

  //! Construct a parse table with the specified keywords.
  Class_Parser_Base(Child &child, Class_Parser_Keyword<Class> const *raw_table,
                    unsigned const count);

  //! This class is meant to be heritable.
  virtual ~Class_Parser_Base(void) = default;

  // MANIPULATORS

  //! Add keywords to the parser.
  void add(Keyword const *table, size_t count) noexcept(false);

  //! Add keywords to the parser. Used for derived parsers.
  template <class Base> void add(std::vector<Class_Parser_Keyword<Base>> const &) noexcept(false);

  //! Remove a keyword from the table.
  void remove(char const *);

  //! Request a change in capacity.
  void reserve(typename std::vector<Keyword>::size_type n) { table_.reserve(n); }

  // ACCESSORS

  //! Return the number of elements in the vector
  // using std::vector<Keyword>::size;
  typename std::vector<Keyword>::size_type size() const { return table_.size(); }

  // SERVICES

  //! Parse a token stream.
  std::shared_ptr<Class> parse(Token_Stream &tokens);

  //! Check the class invariants
  bool check_class_invariants() const;

  // STATIC

  static std::vector<Keyword> const &table() noexcept { return table_; }

protected:
  // IMPLEMENTATION

  virtual void check_completeness(Token_Stream &) = 0;
  virtual Class *create_object() = 0;

private:
  // TYPEDEFS AND ENUMERATIONS

  //----------------------------------------------------------------------------------------------//
  /*!
   * \brief Ordering predicate for Keyword
   *
   * Defines a total ordering for Keyword compatible with STL sort and search routines.
   */

  class Keyword_Compare_ {
  public:
    bool operator()(Keyword const &k1, Keyword const &k2) const;
    bool operator()(Keyword const &keyword, Token const &token) const noexcept;
  };

  // IMPLEMENTATION

  //! Sort and check the table following the addition of new keywords
  void sort_table_() noexcept(false);

  // DATA

  Child &child_;

  // STATIC

  static std::vector<Keyword> table_;
};

//------------------------------------------------------------------------------------------------//
//! Check whether a keyword is well-formed.

template <typename Class> bool Is_Well_Formed_Keyword(Class_Parser_Keyword<Class> const &key);

} // namespace rtt_parser

#include "Class_Parser.i.hh"

#endif // rtt_Class_Parser_HH

//------------------------------------------------------------------------------------------------//
// end of Parse_Table.hh
//------------------------------------------------------------------------------------------------//
