//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   parser/test/tstToken.cc
 * \author Kent G. Budge
 * \date   Feb 18 2003
 * \brief  Unit test for the class rtt_parser::Token.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"
#include "parser/Token.hh"

#ifdef _MSC_VER
#undef ERROR
#endif

using namespace std;
using namespace rtt_parser;
using namespace rtt_dsxx;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

void token_test(UnitTest &ut) {
  Token end_token(END, "parser_test 1");
  Token other_token('$', "parser_test 2");
  Token real_token(REAL, "+1.56e-3", "parser_test 3");

  if (end_token.type() != END)
    ut.failure(__LINE__);
  if (other_token.type() != OTHER)
    ut.failure(__LINE__);
  if (real_token.type() != REAL)
    ut.failure(__LINE__);

  if (other_token.text() != "$")
    ut.failure(__LINE__);
  if (real_token.text() != "+1.56e-3")
    ut.failure(__LINE__);

  if (end_token.location() != "parser_test 1")
    ut.failure(__LINE__);
  if (other_token.location() != "parser_test 2")
    ut.failure(__LINE__);
  if (real_token.location() != "parser_test 3")
    ut.failure(__LINE__);

  if (!Is_Integer_Text("057133"))
    ut.failure(__LINE__);
  if (Is_Integer_Text("08223"))
    ut.failure(__LINE__);
  if (!Is_Integer_Text("663323"))
    ut.failure(__LINE__);
  if (!Is_Integer_Text("0x33a8"))
    ut.failure(__LINE__);
  if (!Is_Keyword_Text("This is a test"))
    ut.failure(__LINE__);
  if (Is_Keyword_Text("This is 1 test"))
    ut.failure(__LINE__);
  if (!Is_Keyword_Text("_"))
    ut.failure(__LINE__);
  if (Is_Keyword_Text("_$"))
    ut.failure(__LINE__);
  if (!Is_Keyword_Text("__"))
    ut.failure(__LINE__);
  if (!Is_Real_Text("+1.56e-3"))
    ut.failure(__LINE__);
#ifndef _MSC_VER
  if (Is_Real_Text("1.39d-3"))
    ut.failure(__LINE__);
#endif
  if (!Is_String_Text("\"This is a test.\""))
    ut.failure(__LINE__);
  if (Is_String_Text("\"This is a test"))
    ut.failure(__LINE__);
  if (Is_String_Text("This is a test"))
    ut.failure(__LINE__);
  if (!Is_String_Text(R"("Backslash \ test")"))
    ut.failure(__LINE__);
  if (Is_String_Text("\"Backslash \\"))
    ut.failure(__LINE__);
  if (Is_Other_Text(""))
    ut.failure(__LINE__);
  if (Is_Other_Text("a"))
    ut.failure(__LINE__);
  if (Is_Other_Text(" "))
    ut.failure(__LINE__);
  if (Is_Other_Text("_"))
    ut.failure(__LINE__);
  if (Is_Other_Text("..."))
    ut.failure(__LINE__);
  if (!Is_Other_Text("=="))
    ut.failure(__LINE__);
  if (!Is_Other_Text("!="))
    ut.failure(__LINE__);
  if (!Is_Other_Text("<="))
    ut.failure(__LINE__);
  if (!Is_Other_Text(">="))
    ut.failure(__LINE__);
  if (!Is_Other_Text("&&"))
    ut.failure(__LINE__);
  if (!Is_Other_Text("||"))
    ut.failure(__LINE__);
  if (Is_Other_Text("!!"))
    ut.failure(__LINE__);

  if (Is_Text_Token(rtt_parser::ERROR))
    ut.failure(__LINE__);
  if (Is_Text_Token(EXIT))
    ut.failure(__LINE__);

  if (Token(REAL, "2", "") == Token(REAL, "3", ""))
    ut.failure(__LINE__);
  if (Token(REAL, "2", "1") == Token(REAL, "2", "3"))
    ut.failure(__LINE__);

  if (real_token == end_token) {
    ut.failure("unlike token equality test did NOT return false");
  } else {
    ut.passes("unlike token equality test returned false");
  }

  // Can only be true? clang-tidy misc-redundant-expression
  //
  // if (real_token == real_token) {
  //   ut.passes("like token equality test did returned true");
  // } else {
  //   ut.failure("like token equality test did NOT return true");
  // }

  return;
}

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
