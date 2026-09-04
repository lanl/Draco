//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   parser/test/tstAbstract_Class_Parser.cc
 * \author Kent G. Budge
 * \brief  Test the Abstract_Class_Parser template
 * \note   Copyright (C) 2010-2026 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "dsxx/Release.hh"
#include "dsxx/ScalarUnitTest.hh"
#include "dsxx/Soft_Equivalence.hh"
#include "dsxx/UnitTest.hh"
#include "parser/Abstract_Class_Parser.hh"
#include "parser/Class_Parse_Table.hh"
#include "parser/Expression.hh"
#include "parser/File_Token_Stream.hh"
#include "parser/Parse_Table.hh"
#include "parser/String_Token_Stream.hh"
#include "parser/Token_Stream.hh"
#include "parser/utilities.hh"
#include <array>
#include <iostream>
#include <limits>
#include <memory>
#include <string>

using namespace std;
using namespace rtt_dsxx;
using namespace rtt_parser;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

// The following would typically be declared in a file associated with the Parent class.
namespace {
class Parent {
public:
  virtual ~Parent() = default;
  Parent() = default;
  Parent(Parent const &rhs) = delete;
  Parent(Parent &&rhs) noexcept = delete;
  Parent &operator=(Parent const &rhs) = delete;
  Parent &operator=(Parent &&rhs) noexcept = delete;
  virtual string name() = 0;
};
} // namespace

//------------------------------------------------------------------------------------------------//
namespace rtt_parser {
template <> class Class_Parse_Table<Parent> {
public:
  // TYPEDEFS

  using Return_Class = Parent;

  // MANAGEMENT

  Class_Parse_Table() {
    child_.reset();
    current_ = this;
  }

  // SERVICES

  [[nodiscard]] Parse_Table const &parse_table() const { return parse_table_; }

  [[nodiscard]] bool allow_exit() const { return false; }

  void check_completeness(rtt_parser::Token_Stream &tokens) {
    if (child_ == std::shared_ptr<Parent>()) {
      tokens.report_semantic_error("no parent specified");
    }
  }

  std::shared_ptr<Parent> create_object() { return child_; }

  // STATICS

  static void register_model(string const &keyword,
                             std::shared_ptr<Parent> parse_function(rtt_parser::Token_Stream &)) {
    Abstract_Class_Parser<Parent, get_parse_table_, get_parsed_object_>::register_child(
        keyword, parse_function);
  }

private:
  // IMPLEMENTATION

  static Parse_Table &get_parse_table_() { return parse_table_; }

  static std::shared_ptr<Parent> &get_parsed_object_() { return child_; }

  // DATA

  static std::shared_ptr<Parent> child_;
  static Class_Parse_Table *current_;
  static Parse_Table parse_table_;
};

std::shared_ptr<Parent> Class_Parse_Table<Parent>::child_;
Class_Parse_Table<Parent> *Class_Parse_Table<Parent>::current_;
Parse_Table Class_Parse_Table<Parent>::parse_table_;

//------------------------------------------------------------------------------------------------//
// Specialization of the parse_class function for T=Parent
template <> std::shared_ptr<Parent> parse_class<Parent>(rtt_parser::Token_Stream &tokens) {
  return parse_class_from_table<Class_Parse_Table<Parent>>(tokens);
}

} // namespace rtt_parser

//------------------------------------------------------------------------------------------------//
// The following is what you would expect to find in a file associated with the Son class.
namespace {
class Son : public Parent {
public:
  string name() override { return "son"; }

  explicit Son(double /*snip_and_snails*/) noexcept {}
};
} // namespace

namespace { // Place this variable in an unnamed namespace instead of using 'static'.
double parsed_snips_and_snails;

void parse_snips_and_snails(rtt_parser::Token_Stream &tokens, int /*unused*/) {
  if (parsed_snips_and_snails >= 0.0) {
    tokens.report_semantic_error("snips and snails already specified");
  }

  parsed_snips_and_snails = parse_real(tokens);
  if (parsed_snips_and_snails < 0.0) {
    tokens.report_semantic_error("snips and snails must not be "
                                 "negative");
    parsed_snips_and_snails = 2;
  }
}

} // namespace

namespace rtt_parser {
template <> class Class_Parse_Table<Son> {
public:
  // TYPEDEFS

  using Return_Class = Son;

  // MANAGEMENT

  Class_Parse_Table() {
    if (!parse_table_is_initialized_) {
      std::array<Keyword, 1> const keywords{
          Keyword{"snips and snails", parse_snips_and_snails, 0, ""}};
      parse_table_.add(keywords.data(), keywords.size());
      parse_table_is_initialized_ = true;
    }

    parsed_snips_and_snails = -1;
    current_ = this;
  }

  // SERVICES

  [[nodiscard]] bool allow_exit() const { return false; }

  [[nodiscard]] Parse_Table const &parse_table() const { return parse_table_; }

  void check_completeness(rtt_parser::Token_Stream &tokens) {
    if (rtt_dsxx::soft_equiv(parsed_snips_and_snails, -1.0,
                             std::numeric_limits<double>::epsilon())) {
      tokens.report_semantic_error("no snips and snails specified");
    }
  }

  std::shared_ptr<Son> create_object() {
    std::shared_ptr<Son> Result = std::make_shared<Son>(parsed_snips_and_snails);
    return Result;
  }

private:
  // STATIC

  static Class_Parse_Table *current_;
  static Parse_Table parse_table_;
  static bool parse_table_is_initialized_;
};

//------------------------------------------------------------------------------------------------//
Class_Parse_Table<Son> *Class_Parse_Table<Son>::current_;
Parse_Table Class_Parse_Table<Son>::parse_table_;
bool Class_Parse_Table<Son>::parse_table_is_initialized_ = false;

//------------------------------------------------------------------------------------------------//
template <> std::shared_ptr<Son> parse_class<Son>(rtt_parser::Token_Stream &tokens) {
  return parse_class_from_table<Class_Parse_Table<Son>>(tokens);
}

} // end namespace rtt_parser

//------------------------------------------------------------------------------------------------//
// The following is what you would expect to find in a file associated with the Daughter class.

namespace {

class Daughter : public Parent {
public:
  string name() override { return "daughter"; }

  explicit Daughter(double /*sugar_and_spice*/) noexcept {}
};

double parsed_sugar_and_spice;

void parse_sugar_and_spice(rtt_parser::Token_Stream &tokens, int /*unused*/) {
  if (parsed_sugar_and_spice >= 0.0) {
    tokens.report_semantic_error("sugar and spice already specified");
  }

  parsed_sugar_and_spice = parse_real(tokens);
  if (parsed_sugar_and_spice < 0.0) {
    tokens.report_semantic_error("sugar and spice must not be negative");
    parsed_sugar_and_spice = 2;
  }
}

} // namespace

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

    parsed_sugar_and_spice = -1;
    current_ = this;
  }

  // SERVICES

  [[nodiscard]] bool allow_exit() const { return false; }

  [[nodiscard]] Parse_Table const &parse_table() const { return parse_table_; }

  void check_completeness(rtt_parser::Token_Stream &tokens) {
    if (rtt_dsxx::soft_equiv(parsed_sugar_and_spice, -1.0,
                             std::numeric_limits<double>::epsilon())) {
      tokens.report_semantic_error("no sugar and spice specified");
    }
  }

  std::shared_ptr<Daughter> create_object() {
    std::shared_ptr<Daughter> Result = std::make_shared<Daughter>(parsed_sugar_and_spice);
    return Result;
  }

private:
  // STATIC

  static Class_Parse_Table *current_;
  static Parse_Table parse_table_;
  static bool parse_table_is_initialized_;
};

Class_Parse_Table<Daughter> *Class_Parse_Table<Daughter>::current_;
Parse_Table Class_Parse_Table<Daughter>::parse_table_;
bool Class_Parse_Table<Daughter>::parse_table_is_initialized_;

//------------------------------------------------------------------------------------------------//
template <> std::shared_ptr<Daughter> parse_class<Daughter>(rtt_parser::Token_Stream &tokens) {
  return parse_class_from_table<Class_Parse_Table<Daughter>>(tokens);
}

} // end namespace rtt_parser

//------------------------------------------------------------------------------------------------//
/* the following would typically live in some client file for Parent. */

namespace { // Place this variable in an unnamed namespace instead of using 'static'.

std::shared_ptr<Parent> parent;

void parse_parent(rtt_parser::Token_Stream &tokens, int /*unused*/) {
  parent = parse_class<Parent>(tokens);
}

std::array<Keyword, 1> const top_keywords{Keyword{"parent", parse_parent, 0, ""}};
Parse_Table top_parse_table(top_keywords.data(), top_keywords.size());

std::shared_ptr<Parent> parse_son(rtt_parser::Token_Stream &tokens) {
  return parse_class<Son>(tokens);
}

std::shared_ptr<Parent> parse_daughter(rtt_parser::Token_Stream &tokens) {
  return parse_class<Daughter>(tokens);
}

//------------------------------------------------------------------------------------------------//
/* Test all the above */

void test(UnitTest &ut) {
  Class_Parse_Table<Parent>::register_model("son", parse_son);

  Class_Parse_Table<Parent>::register_model("daughter", parse_daughter);

  // Build path for the input file
  string const sadInputFile(ut.getTestSourcePath() + std::string("sons_and_daughters.inp"));

  File_Token_Stream tokens(sadInputFile);

  top_parse_table.parse(tokens);

  cout << parent->name() << '\n';

  if (tokens.error_count() == 0 && parent != std::shared_ptr<Parent>() && parent->name() == "son") {
    PASSMSG("Parsed son correctly");
  } else {
    FAILMSG("Did NOT parse son correctly");
  }
}

} // namespace

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  ScalarUnitTest ut(argc, argv, release);
  try {
    test(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstAbstract_Class_Parser.cc
//------------------------------------------------------------------------------------------------//
