//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   parser/test/tstClass_Parser.cc
 * \author Kent Budge
 * \date   Mon Aug 28 07:36:50 2006
 * \note   Copyright (C) 2019-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"
#include "parser/Class_Parse_Table.hh"
#include "parser/String_Token_Stream.hh"
#include "parser/utilities.hh"

using namespace std;
using namespace rtt_dsxx;
using namespace rtt_parser;

//------------------------------------------------------------------------------------------------//
class DummyClass {
public:
  explicit DummyClass(double const insouciance_in) noexcept : insouciance(insouciance_in) {}

  double Get_Insouciance() const { return insouciance; }

private:
  double insouciance;
};

namespace rtt_parser {
//------------------------------------------------------------------------------------------------//
template <> class Class_Parse_Table<DummyClass> : public Class_Parse_Table_Base<DummyClass, false> {
public:
  // MANAGEMENT

  explicit Class_Parse_Table(bool context = false);

  // SERVICES

  void check_completeness(Token_Stream &tokens);

  std::shared_ptr<DummyClass> create_object();

protected:
  // DATA

  double parsed_insouciance{-1.0};
  bool context;

private:
  // IMPLEMENTATION

  static void parse_insouciance_(Token_Stream &tokens, int /*unused*/);
};

//------------------------------------------------------------------------------------------------//
template <>
std::shared_ptr<DummyClass> parse_class<DummyClass>(Token_Stream &tokens, bool const &context) {
  return parse_class_from_table<Class_Parse_Table<DummyClass>>(tokens, context);
}

//------------------------------------------------------------------------------------------------//
template <> std::shared_ptr<DummyClass> parse_class<DummyClass>(Token_Stream &tokens) {
  return parse_class_from_table<Class_Parse_Table<DummyClass>>(tokens, false);
}

//------------------------------------------------------------------------------------------------//
void Class_Parse_Table<DummyClass>::parse_insouciance_(Token_Stream &tokens, int /*unused*/) {
  tokens.check_semantics(current_->parsed_insouciance < 0.0,
                         "duplicate specification of insouciance");

  current_->parsed_insouciance = parse_real(tokens);
  if (current_->parsed_insouciance < 0) {
    tokens.report_semantic_error("insouciance must be nonnegative");
    current_->parsed_insouciance = 1;
  }

  if (current_->context) {
    current_->parsed_insouciance = -current_->parsed_insouciance;
  }
}

//------------------------------------------------------------------------------------------------//
Class_Parse_Table<DummyClass>::Class_Parse_Table(bool const context_in)
    : context(context_in) // sentinel value
{
  std::array<Keyword, 1> const keywords{Keyword{"insouciance", parse_insouciance_, 0, ""}};
  initialize(keywords.data(), sizeof(keywords));
}

//------------------------------------------------------------------------------------------------//
void Class_Parse_Table<DummyClass>::check_completeness(Token_Stream &tokens) {
  tokens.check_semantics(context || parsed_insouciance >= 0, "insouciance was not specified");
}

//------------------------------------------------------------------------------------------------//
std::shared_ptr<DummyClass> Class_Parse_Table<DummyClass>::create_object() {
  std::shared_ptr<DummyClass> Result = std::make_shared<DummyClass>(parsed_insouciance);
  return Result;
}

} // namespace rtt_parser

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

void tstClass_Parser(UnitTest &ut) {
  double const eps = std::numeric_limits<double>::epsilon();
  string text = "insouciance = 3.3\nend\n";
  String_Token_Stream tokens(text);

  std::shared_ptr<DummyClass> dummy = parse_class<DummyClass>(tokens);

  ut.check(dummy != nullptr, "parsed the class object", true);
  ut.check(rtt_dsxx::soft_equiv(dummy->Get_Insouciance(), 3.3, eps),
           "parsed the insouciance correctly");

  tokens.rewind();
  dummy = parse_class<DummyClass>(tokens, true);

  ut.check(dummy != nullptr, "parsed the indolent class object", true);
  ut.check(rtt_dsxx::soft_equiv(dummy->Get_Insouciance(), -3.3, eps),
           "parsed the indolent insouciance correctly");

  // Test that missing end is caught.

  text = "insouciance = 3.3\n";
  String_Token_Stream etokens(text);

  bool good = false;
  try {
    std::shared_ptr<DummyClass> dum = parse_class<DummyClass>(etokens);
  } catch (Syntax_Error &) {
    good = true;
  }
  ut.check(good, "catches missing end");

  tokens.rewind();
  good = false;
  try {
    std::shared_ptr<DummyClass> dum = parse_class<DummyClass>(etokens, true);
  } catch (Syntax_Error &) {
    good = true;
  }
  ut.check(good, "indolent catches missing end");
}

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  ScalarUnitTest ut(argc, argv, release);
  try {
    tstClass_Parser(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstClass_Parser.cc
//------------------------------------------------------------------------------------------------//
