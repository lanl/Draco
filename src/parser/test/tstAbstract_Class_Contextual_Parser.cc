//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   parser/test/tstAbstract_Class_Contextual_Parser.cc
 * \author Kent G. Budge
 * \date   Tue Nov  9 14:34:11 2010
 * \brief  Test the Abstract_Class_Contextual_Parser template
 * \note   Copyright (C) 2016-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"
#include "parser/Abstract_Class_Parser.hh"
#include "parser/Class_Parse_Table.hh"
#include "parser/File_Token_Stream.hh"
#include "parser/utilities.hh"

using namespace std;
using namespace rtt_dsxx;
using namespace rtt_parser;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//
/*
 * Declare an abstract class, Parent, for which we wish to write a constructor. The following would
 * typically be declared in a file named Parent.hh
 */
class Parent {
public:
  explicit Parent(int const magic) : magic_(magic) {}
  // The argument is provided by the parser context and has a magic value, 42, if the parser has
  // handled context correctly.

  int magic() const { return magic_; }

  virtual ~Parent() = default;

  Parent(Parent const &rhs) = delete;
  Parent(Parent &&rhs) noexcept = delete;
  Parent &operator=(Parent const &rhs) = delete;
  Parent &operator=(Parent &&rhs) noexcept = delete;

  virtual string name() = 0;
  // Makes this class abstract, and gives us a way to test later which child has actually been
  // parsed.

private:
  int magic_;
};

//------------------------------------------------------------------------------------------------//
/*
 * Declare a parse table for parsing objects derived from Parent. Note that this must be in the
 * rtt_parser namespace regardless of which namespace Parent lives in.
 */
namespace rtt_parser {

template <> class Class_Parse_Table<Parent> {
public:
  // TYPEDEFS

  using Return_Class = Parent;

  // MANAGEMENT

  explicit Class_Parse_Table(int const &context_in) : context(context_in) {
    child_.reset();
    // Should contain a pointer to a child of Parent after parsing is complete.
    current_ = this;
    // If there are any parsed parameters common to all children of Parent, they should be placed in
    // Class_Parse_Table<Parent> along with the associated parsing routines, which will populate the
    // parameters through the current_ pointer. This is necessary becuase the parsing routines must
    // be static functions. (At least until we get ambitious with closures.
  }

  // SERVICES

  bool allow_exit() const { return false; }

  void check_completeness(Token_Stream &tokens) {
    if (!child_) {
      tokens.report_semantic_error("no parent specified");
    }
  }

  std::shared_ptr<Parent> create_object() { return child_; }

  // STATICS

  Parse_Table const &parse_table() const { return parse_table_; }

  static void register_model(string const &keyword,
                             std::shared_ptr<Parent> parse_function(Token_Stream &, int const &)) {
    Abstract_Class_Parser<Parent, get_parse_table_, get_parsed_object_,
                          Contextual_Parse_Functor<Parent, int, get_context_>>::
        register_child(keyword,
                       Contextual_Parse_Functor<Parent, int, get_context_>(parse_function));
    // This function allows downstream developers to add new daughter classes to the parser without
    // having to modify upstream code.
  }

protected:
  // DATA

  // This is where any parameters shared by all children of Parent would be placed.

  int context;
  // This is the parser context.

private:
  // IMPLEMENTATION

  static Parse_Table &get_parse_table_() { return parse_table_; }

  static std::shared_ptr<Parent> &get_parsed_object_() { return child_; }

  static int const &get_context_() { return current_->context; }

  // DATA

  static std::shared_ptr<Parent> child_;
  static Class_Parse_Table *current_;
  static Parse_Table parse_table_;
};

std::shared_ptr<Parent> Class_Parse_Table<Parent>::child_;
Class_Parse_Table<Parent> *Class_Parse_Table<Parent>::current_;
Parse_Table Class_Parse_Table<Parent>::parse_table_;

//------------------------------------------------------------------------------------------------//
// Specialization of the parse_class template function for T=Parent
template <> std::shared_ptr<Parent> parse_class(Token_Stream &tokens, int const &context) {
  return parse_class_from_table<Class_Parse_Table<Parent>>(tokens, context);
}

} // namespace rtt_parser

//------------------------------------------------------------------------------------------------//
/*
 * We now declare a child class derived from Parent.  The following would normally be placed in the
 * file Son.hh
 */
class Son : public Parent {
public:
  string name() override { return "son"; }

  Son(double /*snip_and_snails*/, int const context) : Parent(context) {}
  // "snips and snails" is provided by the parser based on the parsed specification, while context
  // is provided by the parser context provided by the client. Since this is a toy example, we don't
  // actually use the parameter, but we do want to check that the context is got right.
};

//------------------------------------------------------------------------------------------------//
/*
 * Now declare a parser class for parsing specifications for Son. Typically the child parser class
 * will be derived from the parent parser class so that parse code for common parameters does not
 * have to be duplicated.
 */
namespace rtt_parser {
template <> class Class_Parse_Table<Son> : public Class_Parse_Table<Parent> {
public:
  // TYPEDEFS

  using Return_Class = Son;

  // MANAGEMENT

  explicit Class_Parse_Table(int const context_in) : Class_Parse_Table<Parent>(context_in) {
    if (!parse_table_is_initialized_) {
      // The parser class must populate the parse_table_ with the keywords and parse functions
      // needed to parse a specification. This is done once the first time any
      // Class_Parse_Table<Son> object is constructed.

      std::array<Keyword, 1> const keywords{
          Keyword{"snips and snails", parse_snips_and_snails, 0, ""}};
      parse_table_.add(keywords.data(), keywords.size());
      parse_table_is_initialized_ = true;
    }

    // Initialize the parameters about to be parsed, typical with sentinel values that indicate
    // whether the parameter has been found in the parsed specification.
    snips_and_snails = -1;

    // Set the current_ pointer to this object so that the parse routines, which must be static,
    // will know where to find the Class_Parse_Table<Son> object in which to store parsed
    // parameters.
    current_ = this;
    // If the parser is meant to support reentrancy, then the old current_ pointer needs to be saved
    // somewhere (such as a stack) so it can be retrieved later.
  }

  // SERVICES

  bool allow_exit() const { return false; }

  Parse_Table const &parse_table() const { return parse_table_; }

  void check_completeness(Token_Stream &tokens) {
    if (rtt_dsxx::soft_equiv(snips_and_snails, -1.0, std::numeric_limits<double>::epsilon())) {
      tokens.report_semantic_error("no snips and snails specified");
    }
  }

  std::shared_ptr<Son> create_object() {
    std::shared_ptr<Son> Result(new Son(snips_and_snails, context));
    return Result;
  }

protected:
  // Parameter to be parsed

  double snips_and_snails;

private:
  // STATIC

  static void parse_snips_and_snails(Token_Stream &tokens, int /*unused*/) {
    if (current_->snips_and_snails >= 0.0) {
      tokens.report_semantic_error("snips and snails already specified");
    }

    current_->snips_and_snails = parse_real(tokens);
    if (current_->snips_and_snails < 0.0) {
      tokens.report_semantic_error("snips and snails must not be "
                                   "negative");
      current_->snips_and_snails = 2;
      // It's customary to set parameters that have an invalid value specified to some benign valid
      // value.
    }
  }

  static Class_Parse_Table *current_;
  static Parse_Table parse_table_;
  static bool parse_table_is_initialized_;
};

//------------------------------------------------------------------------------------------------//
Class_Parse_Table<Son> *Class_Parse_Table<Son>::current_;
Parse_Table Class_Parse_Table<Son>::parse_table_;
bool Class_Parse_Table<Son>::parse_table_is_initialized_ = false;

//------------------------------------------------------------------------------------------------//
template <> std::shared_ptr<Son> parse_class<Son>(Token_Stream &tokens, int const &context) {
  return parse_class_from_table<Class_Parse_Table<Son>>(tokens, context);
}

} // end namespace rtt_parser

//------------------------------------------------------------------------------------------------//
/*
 * Now define a second child of Parent, which we will (whimsically) call Daughter. The following
 * would typicall be placed in the file Daughter.hh
 */
class Daughter : public Parent {
public:
  string name() override { return "daughter"; }

  explicit Daughter(double /*sugar_and_spice*/) : Parent(0) {}
  // This child doesn't care about the context, which is perfectly acceptable (if it makes sense).
};

//------------------------------------------------------------------------------------------------//
/*
 * Define a parser class now for Daughter. This is similar to what we do for Son so we will go light
 * on comments.
 */
namespace rtt_parser {
template <> class Class_Parse_Table<Daughter> {
public:
  // TYPEDEFS

  using Return_Class = Daughter;

  // MANAGEMENT

  Class_Parse_Table() {
    if (!parse_table_is_initialized_) {
      std::array<Keyword, 1> const keywords{
          Keyword{"sugar and spice", parse_sugar_and_spice, 0, ""}};
      parse_table_.add(keywords.data(), keywords.size());
      parse_table_is_initialized_ = true;
    }

    sugar_and_spice = -1;

    current_ = this;
  }

  // SERVICES

  bool allow_exit() const { return false; }

  Parse_Table const &parse_table() const { return parse_table_; }

  void check_completeness(Token_Stream &tokens) {
    if (rtt_dsxx::soft_equiv(sugar_and_spice, -1.0, std::numeric_limits<double>::epsilon())) {
      tokens.report_semantic_error("no sugar and spice specified");
    }
  }

  std::shared_ptr<Daughter> create_object() {
    std::shared_ptr<Daughter> Result(new Daughter(sugar_and_spice));
    return Result;
  }

protected:
  double sugar_and_spice;

private:
  // STATIC

  static void parse_sugar_and_spice(Token_Stream &tokens, int /*unused*/) {
    if (current_->sugar_and_spice >= 0.0) {
      tokens.report_semantic_error("sugar and spice already specified");
    }

    current_->sugar_and_spice = parse_real(tokens);
    if (current_->sugar_and_spice < 0.0) {
      tokens.report_semantic_error("sugar and spice must not be "
                                   "negative");
      current_->sugar_and_spice = 2;
    }
  }

  static Class_Parse_Table *current_;
  static Parse_Table parse_table_;
  static bool parse_table_is_initialized_;
};

Class_Parse_Table<Daughter> *Class_Parse_Table<Daughter>::current_;
Parse_Table Class_Parse_Table<Daughter>::parse_table_;
bool Class_Parse_Table<Daughter>::parse_table_is_initialized_;

//------------------------------------------------------------------------------------------------//
template <> std::shared_ptr<Daughter> parse_class<Daughter>(Token_Stream &tokens) {
  return parse_class_from_table<Class_Parse_Table<Daughter>>(tokens);
}

} // end namespace rtt_parser

// The following are the kinds of parse functions that a downwind developer might write for his own
// children of the Parent class.

std::shared_ptr<Parent> parse_son(Token_Stream &tokens, int const &context) {
  return parse_class<Son>(tokens, context);
}

std::shared_ptr<Parent> parse_daughter(Token_Stream &tokens, int const & /*unused*/) {
  return parse_class<Daughter>(tokens);
}

//------------------------------------------------------------------------------------------------//
/* Test all the above */

void test(UnitTest &ut) {
  // Register the children parse functions with the Parent parse class.
  Class_Parse_Table<Parent>::register_model("son", parse_son);
  Class_Parse_Table<Parent>::register_model("daughter", parse_daughter);

  // Build path for the input file containing a test specification for a Son.
  string const sadInputFile(ut.getTestSourcePath() +
                            std::string("contextual_sons_and_daughters.inp"));

  File_Token_Stream tokens(sadInputFile);

  std::shared_ptr<Parent> parent = parse_class<Parent>(tokens, 42);
  // We choose 42 as the magic value for our context, in honor of hitchikers across the Galaxy.

  cout << parent->name() << endl;

  ut.check(tokens.error_count() == 0, "parsed without error");
  ut.check(parent != nullptr, "created parent", true);
  ut.check(parent->name() == "son", "parent is son");
  ut.check(parent->magic() == 42, "context");
}

//------------------------------------------------------------------------------------------------//

int main(int argc, char *argv[]) {
  ScalarUnitTest ut(argc, argv, release);
  try {
    test(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstAbstract_Class_Contextual_Parser.cc
//------------------------------------------------------------------------------------------------//
