//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   parser/test/tstToken.cc
 * \author Kent G. Budge
 * \brief  Unit test for the class rtt_parser::Token.
 * \note   Copyright (C) 2010-2026 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "dsxx/Release.hh"
#include "dsxx/ScalarUnitTest.hh"
#include "dsxx/UnitTest.hh"
#include "parser/Token.hh"

#ifdef _MSC_VER
#undef ERROR
#endif

namespace {

using namespace std;
using namespace rtt_parser;
using namespace rtt_dsxx;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

void token_test(UnitTest &ut) {
  Token const end_token(END, "parser_test 1");
  Token const other_token('$', "parser_test 2");
  Token const real_token(REAL, "+1.56e-3", "parser_test 3");

  if (end_token.type() != END)
    FAILMSG(__LINE__);
  if (other_token.type() != OTHER)
    FAILMSG(__LINE__);
  if (real_token.type() != REAL)
    FAILMSG(__LINE__);

  if (other_token.text() != "$")
    FAILMSG(__LINE__);
  if (real_token.text() != "+1.56e-3")
    FAILMSG(__LINE__);

  if (end_token.location() != "parser_test 1")
    FAILMSG(__LINE__);
  if (other_token.location() != "parser_test 2")
    FAILMSG(__LINE__);
  if (real_token.location() != "parser_test 3")
    FAILMSG(__LINE__);

  if (!Is_Integer_Text("057133"))
    FAILMSG(__LINE__);
  if (Is_Integer_Text("08223"))
    FAILMSG(__LINE__);
  if (!Is_Integer_Text("663323"))
    FAILMSG(__LINE__);
  if (!Is_Integer_Text("0x33a8"))
    FAILMSG(__LINE__);
  if (!Is_Keyword_Text("This is a test"))
    FAILMSG(__LINE__);
  if (Is_Keyword_Text("This is 1 test"))
    FAILMSG(__LINE__);
  if (!Is_Keyword_Text("_"))
    FAILMSG(__LINE__);
  if (Is_Keyword_Text("_$"))
    FAILMSG(__LINE__);
  if (!Is_Keyword_Text("__"))
    FAILMSG(__LINE__);
  if (!Is_Real_Text("+1.56e-3"))
    FAILMSG(__LINE__);
#ifndef _MSC_VER
  if (Is_Real_Text("1.39d-3"))
    FAILMSG(__LINE__);
#endif
  if (!Is_String_Text("\"This is a test.\""))
    FAILMSG(__LINE__);
  if (Is_String_Text("\"This is a test"))
    FAILMSG(__LINE__);
  if (Is_String_Text("This is a test"))
    FAILMSG(__LINE__);
  if (!Is_String_Text(R"("Backslash \ test")"))
    FAILMSG(__LINE__);
  if (Is_String_Text("\"Backslash \\"))
    FAILMSG(__LINE__);
  if (Is_Other_Text(""))
    FAILMSG(__LINE__);
  if (Is_Other_Text("a"))
    FAILMSG(__LINE__);
  if (Is_Other_Text(" "))
    FAILMSG(__LINE__);
  if (Is_Other_Text("_"))
    FAILMSG(__LINE__);
  if (Is_Other_Text("..."))
    FAILMSG(__LINE__);
  if (!Is_Other_Text("=="))
    FAILMSG(__LINE__);
  if (!Is_Other_Text("!="))
    FAILMSG(__LINE__);
  if (!Is_Other_Text("<="))
    FAILMSG(__LINE__);
  if (!Is_Other_Text(">="))
    FAILMSG(__LINE__);
  if (!Is_Other_Text("&&"))
    FAILMSG(__LINE__);
  if (!Is_Other_Text("||"))
    FAILMSG(__LINE__);
  if (Is_Other_Text("!!"))
    FAILMSG(__LINE__);

  if (Is_Text_Token(rtt_parser::ERROR))
    FAILMSG(__LINE__);
  if (Is_Text_Token(EXIT))
    FAILMSG(__LINE__);

  if (Token(REAL, "2", "") == Token(REAL, "3", ""))
    FAILMSG(__LINE__);
  if (Token(REAL, "2", "1") == Token(REAL, "2", "3"))
    FAILMSG(__LINE__);

  if (real_token == end_token) {
    FAILMSG("unlike token equality test did NOT return false");
  } else {
    PASSMSG("unlike token equality test returned false");
  }

  return;
}

} // namespace

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  ScalarUnitTest ut(argc, argv, release);
  try {
    token_test(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstToken.cc
//------------------------------------------------------------------------------------------------//
