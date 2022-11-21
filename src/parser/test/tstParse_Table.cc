//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   parser/test/tstParse_Table.cc
 * \author Kent G. Budge
 * \date   Feb 18 2003
 * \brief  Unit tests for the Parse_Table class.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"
#include "parser/File_Token_Stream.hh"
#include "parser/Parse_Table.hh"
#include "parser/String_Token_Stream.hh"
#include <array>
#include <cstring>

#ifdef _MSC_VER
#undef ERROR
#endif

using namespace std;
using namespace rtt_parser;
using namespace rtt_dsxx;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

static std::array<const char *, 3> const color{"BLACK", "BLUE", "BLUE GREEN"};
static std::array<bool, 3> color_set;

static void Parse_Color(Token_Stream & /*unused*/, int i) {
  cout << "You have requested " << color[i] << endl;
  color_set[i] = true;
}

static void Parse_Any_Color(Token_Stream &tokens, int /*unused*/) {
  Token token = tokens.shift();
  for (unsigned i = 0; i < sizeof(color) / sizeof(const char *); i++)
    if (!strcmp(token.text().c_str(), color[i])) {
      cout << "You have requested " << color[i] << endl;
      color_set[i] = true;
      return;
    }

  tokens.report_syntax_error(token, "expected a color");
}

std::array<Keyword, 6> const raw_table{
    Keyword{"BLUE", Parse_Color, 1, "main"},
    Keyword{"BLACK", Parse_Color, 0, "main"},
    Keyword{"BLUE GREEN", Parse_Color, 2, "main"},
    Keyword{"BLUISH GREEN", Parse_Color, 2, "main"},
    Keyword{"lower blue", Parse_Color, 2, "main", "keyword to test case sensitivity"},
    Keyword{"COLOR", Parse_Any_Color, 0, "main"}};

std::array<Keyword, 2> const raw_table_2{Keyword{"BLUE", Parse_Color, 1, "main"},
                                         Keyword{"BLACK", Parse_Color, 0, "main"}};

class Error_Token_Stream : public Token_Stream {
public:
  void rewind() override {}

  void comment(string const & /*err*/) override {
    cout << "comment reported to Error_Token_Stream" << endl;
  }

protected:
  void report(Token const & /*token*/, string const & /*err*/) override {
    cout << "error reported to Error_Token_Stream" << endl;
  }

  void report(string const & /*err*/) override {
    cout << "error reported to Error_Token_Stream" << endl;
  }

  Token fill_() override { return {rtt_parser::ERROR, "error"}; }
};

class Colon_Token_Stream : public Token_Stream {
public:
  Colon_Token_Stream() = default;

  void rewind() override {}

  void comment(string const & /*err*/) override {
    cout << "comment reported to Colon_Token_Stream" << endl;
  }

protected:
  void report(Token const & /*token*/, string const & /*err*/) override {
    cout << "error reported to Colon_Token_Stream" << endl;
  }

  void report(string const & /*err*/) override {
    cout << "error reported to Colon_Token_Stream" << endl;
  }

  Token fill_() override {
    switch (count_++) {
    case 0:
      return {';', ""};
    case 1:
      return {END, "end"};
    case 2:
      return {EXIT, ""};
    default:
      Insist(false, "bad case");
      return {rtt_parser::ERROR, ""}; // dummy return to eliminate warning
    }
  }

private:
  unsigned count_{0};
};

//------------------------------------------------------------------------------------------------//
void tstKeyword(UnitTest &ut) {
  Keyword black = {"BLACK", Parse_Color, 0, "main"};

  if (black == raw_table_2[1])
    PASSMSG("keyword equality operator okay");
  else
    FAILMSG("keyword equality operator NOT okay");

  {
    Keyword key1 = {"BLUE", Parse_Color, 1, "main"};
    Keyword key2 = {"RED", Parse_Color, 1, "main"};
    if (key1 == key2)
      FAILMSG("comparison of dislike keywords FAILED");
  }
  {
    Keyword key1 = {"BLUE", Parse_Color, 1, "main"};
    Keyword key2 = {"BLUE", Parse_Any_Color, 1, "main"};
    if (key1 == key2)
      FAILMSG("comparison of dislike functions FAILED");
  }
  {
    Keyword key1 = {"BLUE", Parse_Color, 1, "main"};
    Keyword key2 = {"BLUE", Parse_Color, 12, "main"};
    if (key1 == key2)
      FAILMSG("comparison of dislike index FAILED");
  }
  {
    Keyword key1 = {"BLUE", Parse_Color, 1, "main"};
    Keyword key2 = {"BLUE", Parse_Color, 1, "second"};
    if (key1 == key2)
      FAILMSG("comparison of dislike module FAILED");
  }

  {
    Keyword key = {nullptr, Parse_Color, 1, "main"};
    if (Is_Well_Formed_Keyword(key))
      FAILMSG("null keyword detect FAILED");
  }
  {
    Keyword key = {"BLUE", nullptr, 1, "main"};
    if (Is_Well_Formed_Keyword(key))
      FAILMSG("null func detect FAILED");
  }
  {
    Keyword key = {".BLUE", Parse_Color, 1, "main"};
    if (Is_Well_Formed_Keyword(key))
      FAILMSG("bad moniker detect FAILED");
  }
  {
    Keyword key = {"_BLUE", Parse_Color, 1, "main"};
    if (!Is_Well_Formed_Keyword(key))
      FAILMSG("moniker with leading underscore FAILED");
  }
  {
    Keyword key = {"BLUE.", Parse_Color, 1, "main"};
    if (Is_Well_Formed_Keyword(key))
      FAILMSG("bad moniker detect FAILED");
  }
}

//------------------------------------------------------------------------------------------------//
void tstParse_Table(UnitTest &ut) {
  Parse_Table table;

  table.reserve(raw_table.size());
  table.add(raw_table.data(), raw_table.size());

  UT_MSG(table.size() == raw_table.size(), "Found expected table size");

  // Build path for the input file "parser_test.inp"
  string const ptInputFile(ut.getTestSourcePath() + std::string("parser_test.inp"));

  File_Token_Stream token_stream(ptInputFile);

  table.parse(token_stream);

  FAIL_IF(!color_set[1]);
  FAIL_IF(token_stream.error_count() != 5);

  token_stream.rewind();

  table.set_flags(Parse_Table::CASE_INSENSITIVE);

  color_set[0] = color_set[1] = false;
  table.parse(token_stream);

  FAIL_IF(!color_set[1]);
  FAIL_IF(token_stream.error_count() != 4);

  {
    String_Token_Stream tokens("BLUE green");
    table.parse(tokens);
    if (tokens.error_count() != 0)
      FAILMSG("Did NOT match mismatched case");
  }
  {
    String_Token_Stream tokens("lower blue");
    table.parse(tokens);
    if (tokens.error_count() != 0)
      FAILMSG("Did NOT match lower case");
  }
  {
    String_Token_Stream tokens("lowe");
    table.parse(tokens);
    if (tokens.error_count() != 1)
      FAILMSG("Did NOT detect partial match case");
  }
  {
    String_Token_Stream tokens("lower bluer");
    table.parse(tokens);
    if (tokens.error_count() != 1)
      FAILMSG("Did NOT detect partial match case");
  }

  token_stream.rewind();

  table.set_flags(Parse_Table::CASE_INSENSITIVE | Parse_Table::PARTIAL_IDENTIFIER_MATCH);

  color_set[0] = color_set[1] = false;
  table.parse(token_stream);

  FAIL_IF(!color_set[1]);
  FAIL_IF(token_stream.error_count() != 3);

  // Test the Get_Flags() function, even if this test is built with the flag
  // --with-dbc=0.
  FAIL_IF(table.get_flags() != 3);

  // Test the check_class_invariants() function, even if this test is built
  // with the flag --with-dbc=0.
  FAIL_IF(!table.check_class_invariants());

  // Check variations on partial match
  {
    String_Token_Stream tokens("BLUEE");
    table.parse(tokens);
    if (tokens.error_count() != 0)
      FAILMSG("Did NOT match partial keyword");
  }
  {
    String_Token_Stream tokens("blue");
    table.parse(tokens);
    if (tokens.error_count() != 0)
      FAILMSG("Did NOT match keyword with wrong case");
  }
  {
    String_Token_Stream tokens("end");
    if (table.parse(tokens).type() != END) {
      FAILMSG("END detection FAILED");
    }
  }
  // Test recovery
  {}

  table.set_flags(Parse_Table::PARTIAL_IDENTIFIER_MATCH);
  {
    String_Token_Stream tokens("BLUEE");
    table.parse(tokens);
    if (tokens.error_count() != 0)
      FAILMSG("Did NOT match partial keyword");
  }
  {
    String_Token_Stream tokens("BLU green");
    table.parse(tokens);
    if (tokens.error_count() != 1)
      FAILMSG("Did NOT detect mismatched case");
  }
  {
    String_Token_Stream tokens("blue");
    table.parse(tokens);
    if (tokens.error_count() != 1)
      FAILMSG("Did NOT detect mismatched case");
  }
  {
    String_Token_Stream tokens("blue green red");
    table.parse(tokens);
    if (tokens.error_count() != 1)
      FAILMSG("Did NOT detect mismatched case");
  }
  {
    String_Token_Stream tokens("BLUE RED");
    table.parse(tokens);
    if (tokens.error_count() != 1)
      FAILMSG("Did NOT detect unknown keyword");
  }
  {
    String_Token_Stream tokens("BLUISH");
    table.parse(tokens);
    if (tokens.error_count() != 1)
      FAILMSG("Did NOT catch partial mismatch");
  }
  {
    String_Token_Stream tokens("end");
    if (table.parse(tokens).type() != END) {
      FAILMSG("END detection FAILED");
    }
    if (table.parse(tokens).type() != EXIT) {
      FAILMSG("exit detection FAILED");
    }
  }
  {
    Colon_Token_Stream tokens;
    tokens.comment("dummy test");
    if (table.parse(tokens).type() != END) {
      FAILMSG("END detection FAILED");
    }
    if (table.parse(tokens).type() != EXIT) {
      FAILMSG("exit detection FAILED");
    }
  }
  // Error handling
  {
    Error_Token_Stream tokens;
    tokens.comment("dummy comment");
    if (table.parse(tokens).type() != rtt_parser::ERROR) {
      FAILMSG("error detection FAILED");
    }
  }

  Parse_Table table_2(raw_table.data(), raw_table.size());

  if (table_2.size() != raw_table.size())
    FAILMSG("test FAILS");

  token_stream.rewind();

  table_2.parse(token_stream);

  if (!color_set[1])
    FAILMSG("test FAILS");

  if (token_stream.error_count() != 5)
    FAILMSG("error count FAILS");

  Keyword test_key = {"THIS SHOULD WORK", Parse_Color, 0, nullptr};
  if (!Is_Well_Formed_Keyword(test_key))
    FAILMSG("test FAILS");

  std::array<Keyword, 2> benign_ambiguous_table{Keyword{"KEY", Parse_Color, 0, nullptr},
                                                Keyword{"KEY", Parse_Color, 0, nullptr}};
  table_2.add(benign_ambiguous_table.data(), benign_ambiguous_table.size());
  token_stream.rewind();
  table_2.parse(token_stream);

  std::array<Keyword, 1> malign_ambiguous_table{Keyword{"KEY", Parse_Color, 1, nullptr}};
  try {
    table_2.add(malign_ambiguous_table.data(), malign_ambiguous_table.size());
    token_stream.rewind();
    table_2.parse(token_stream);
    FAILMSG("did NOT catch ambiguous keyword");
  } catch (invalid_argument const &msg) {
    cout << msg.what() << endl;
    PASSMSG("successfully detected ambiguous keyword");
  }

  // Build path for the input file "recovery.inp"
  string const recInputFile(ut.getTestSourcePath() + std::string("recovery.inp"));

  File_Token_Stream recover_stream(recInputFile);
  table.parse(recover_stream);
  if (recover_stream.error_count() != 2)
    FAILMSG("test FAILS");

  Parse_Table table_3;
  std::array<Keyword, 2> case_ambiguous_table{Keyword{"key", Parse_Color, 0, nullptr},
                                              Keyword{"Key", Parse_Color, 1, nullptr}};
  try {
    table_3.add(case_ambiguous_table.data(), case_ambiguous_table.size());
    table_3.parse(token_stream);
    table_3.set_flags(Parse_Table::CASE_INSENSITIVE);
    token_stream.rewind();
    table_3.parse(token_stream);
    FAILMSG("did NOT catch case-dependent ambiguous keyword");
  } catch (invalid_argument const &msg) {
    cout << msg.what() << endl;
    PASSMSG("successfully detected case-dependent ambiguous keyword");
  }

  Parse_Table table_3a;
  std::array<Keyword, 2> casea_ambiguous_table{Keyword{"key", Parse_Color, 0, nullptr},
                                               Keyword{"Key", Parse_Any_Color, 0, nullptr}};
  try {
    table_3a.add(casea_ambiguous_table.data(), casea_ambiguous_table.size());
    table_3a.parse(token_stream);
    table_3a.set_flags(Parse_Table::CASE_INSENSITIVE);
    token_stream.rewind();
    table_3a.parse(token_stream);
    FAILMSG("did NOT catch case-dependent ambiguous keyword");
  } catch (invalid_argument const &msg) {
    cout << msg.what() << endl;
    PASSMSG("successfully detected case-dependent ambiguous keyword");
  }

  Parse_Table table_4;
  table_4.add(table_4);
  table_4.add(table);

  recover_stream.rewind();
  table_4.parse(recover_stream);
  if (recover_stream.error_count() != 2)
    FAILMSG("test FAILS");

  table.set_flags(Parse_Table::ONCE);
  {
    String_Token_Stream tokens("BLUE, ERROR");
    if (table.parse(tokens).type() != END || tokens.error_count() > 0) {
      FAILMSG("FAILED to end on one token in ONCE mode");
    }
  }

  {
    // Additional test mandated by bug discovery.

    Parse_Table ptable;

    ptable.reserve(raw_table_2.size());
    ptable.add(raw_table_2.data(), raw_table_2.size());

    if (ptable.size() != raw_table_2.size())
      FAILMSG("test FAILS");

    File_Token_Stream ltoken_stream(ptInputFile);

    ptable.parse(ltoken_stream);

    if (ltoken_stream.error_count() != 5)
      FAILMSG("test FAILS");
  }

  // Test keyword removal

  table.remove("BLUE");
  {
    String_Token_Stream tokens("BLUE");
    table.parse(tokens);
    if (tokens.error_count() == 0) {
      FAILMSG("FAILED to remove token");
    }
  }

  return;
}

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  ScalarUnitTest ut(argc, argv, release);
  try {
    tstKeyword(ut);
    tstParse_Table(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstParse_Table.cc
//------------------------------------------------------------------------------------------------//
