//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   parser/test/tstConsole_Token_Stream.cc
 * \author Kent G. Budge
 * \date   Wed May 19 11:26:15 MDT 2004
 * \brief  Unit tests for Console_Token_Stream class.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"
#include "parser/Console_Token_Stream.hh"
#include <sstream>

using namespace std;
using namespace rtt_parser;
using namespace rtt_dsxx;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

void tstConsole_Token_Stream(rtt_dsxx::UnitTest &ut) {
  {
    set<char> ws;
    ws.insert(' ');
    Console_Token_Stream tokens(ws);
    if (tokens.whitespace() != ws)
      FAILMSG("whitespace characters are NOT correctly set");
    else
      PASSMSG("whitespace characters are correctly set");
  }

  {
    Console_Token_Stream tokens;
    tokens.comment("begin Console_Token_Stream tests");
    if (tokens.whitespace() != Text_Token_Stream::default_whitespace)
      FAILMSG("whitespace characters are NOT correct defaults");
    else
      PASSMSG("whitespace characters are correct defaults");

    Token token = tokens.lookahead(3);
    if (token.type() != KEYWORD || token.text() != "COLOR")
      FAILMSG("lookahead(3) does NOT have correct value");
    else
      PASSMSG("lookahead(3) has correct value");

    tokens.report_semantic_error(token, "dummy error");
    if (tokens.error_count() != 1)
      FAILMSG("Dummy error NOT counted properly");
    else
      PASSMSG("Dummy error counted properly");

    tokens.report_semantic_error("dummy error");
    if (tokens.error_count() != 2)
      FAILMSG("Dummy error NOT counted properly");
    else
      PASSMSG("Dummy error counted properly");

    token = tokens.shift();
    if (token.type() != KEYWORD || token.text() != "BLUE")
      FAILMSG("First shift does NOT have correct value");
    else
      PASSMSG("First shift has correct value");

    token = tokens.lookahead();
    if (token.type() != KEYWORD || token.text() != "GENERATE ERROR")
      FAILMSG("Lookahed after first shift does NOT have correct value");
    else
      PASSMSG("Lookahead after first shift has correct value");

    token = tokens.shift();
    if (token.type() != KEYWORD || token.text() != "GENERATE ERROR")
      FAILMSG("Second shift does NOT have correct value");
    else
      PASSMSG("Second shift has correct value");

    token = tokens.shift();
    if (token.type() != KEYWORD || token.text() != "GENERATE ANOTHER ERROR")
      FAILMSG("Third shift does NOT have correct value");
    else
      PASSMSG("Third shift has correct value");

    token = Token('$', "test_parser");
    tokens.pushback(token);

    token = tokens.shift();
    if (token.type() != OTHER || token.text() != "$")
      FAILMSG("Shift after pushback does NOT have correct value");
    else
      PASSMSG("Shift after pushback has correct value");

    bool caught(false);
    try {
      tokens.report_syntax_error(token, "dummy syntax error");
    } catch (const Syntax_Error & /*msg*/) {
      caught = true;
      PASSMSG("Syntax error correctly thrown and caught");
    }
    FAIL_IF_NOT(caught); // FAILMSG("Syntax error NOT correctly thrown");

    token = tokens.shift();
    if (token.type() != KEYWORD || token.text() != "COLOR")
      ITFAILS;

    token = tokens.shift();
    if (token.type() != KEYWORD || token.text() != "BLACK")
      ITFAILS;

    token = tokens.shift();
    if (token.type() != END)
      ITFAILS;

    token = tokens.shift();
    if (token.type() != OTHER || token.text() != "-")
      ITFAILS;

    token = tokens.shift();
    if (token.type() != REAL || token.text() != "1.563e+3")
      ITFAILS;

    token = tokens.shift();
    if (token.type() != REAL || token.text() != "1.563e+3")
      ITFAILS;

    token = tokens.shift();
    if (token.type() != REAL || token.text() != ".563e+3")
      ITFAILS;

    token = tokens.shift();
    if (token.type() != OTHER || token.text() != ".")
      ITFAILS;

    token = tokens.shift();
    if (token.type() != OTHER || token.text() != "-")
      ITFAILS;

    token = tokens.shift();
    if (token.type() != REAL || token.text() != "1.")
      ITFAILS;

    token = tokens.shift();
    if (token.type() != REAL || token.text() != "1.563")
      ITFAILS;

    token = tokens.shift();
    if (token.type() != REAL || token.text() != "1.e+3")
      ITFAILS;

    token = tokens.shift();
    if (token.type() != REAL || token.text() != "1.e3")
      ITFAILS;

    token = tokens.shift();
    if (token.type() != REAL || token.text() != "1e+3")
      ITFAILS;

    token = tokens.shift();
    if (token.type() != REAL || token.text() != "1e3")
      ITFAILS;

    token = tokens.shift();
    if (token.type() != INTEGER || token.text() != "19090")
      ITFAILS;

    token = tokens.shift();
    if (token.type() != INTEGER || token.text() != "01723")
      ITFAILS;

    token = tokens.shift();
    if (token.type() != INTEGER || token.text() != "0x1111a")
      ITFAILS;

    token = tokens.shift();
    if (token.type() != INTEGER || token.text() != "0")
      ITFAILS;

    token = tokens.shift();
    if (token.type() != INTEGER || token.text() != "8123")
      ITFAILS;

    token = tokens.shift();
    if (token.type() != STRING || token.text() != "\"manifest string\"")
      ITFAILS;

    token = tokens.shift();
    FAIL_IF(token.type() != STRING || token.text() != R"("manifest \"string\"")");

    token = tokens.shift();
    if (token.type() != OTHER || token.text() != "@")
      ITFAILS;

    token = tokens.shift();
    if (token.type() != INTEGER || token.text() != "1")
      ITFAILS;

    token = tokens.shift();
    if (token.type() != KEYWORD || token.text() != "e")
      ITFAILS;

    token = tokens.shift();
    if (token.type() != INTEGER || token.text() != "0")
      ITFAILS;

    token = tokens.shift();
    if (token.type() != KEYWORD || token.text() != "x")
      ITFAILS;

    // Test that #include is treated as a syntax error (not supported)
    try {
      tokens.shift();
      ut.failure("Did NOT correctly report #include as error");
    } catch (const Syntax_Error & /*msg*/) {
      PASSMSG("#include not supported error correctly thrown and caught");
    }

    token = tokens.shift();
    if (token.type() != EXIT)
      ITFAILS;
    token = tokens.shift();
    if (token.type() != EXIT)
      ITFAILS;

    // [2015-09-01 KT] After talking to Kent about this implementation, we
    // decided that it does not make sense to rewind an interactive standard
    // input buffer.  In fact, if this is done under an MPI environment
    // (e.g. mpirun -np 1, aprun -n 1, etc. ), the seekg() will return an
    // error condition.

    // {
    //     tokens.rewind();
    //     token = tokens.lookahead();
    //     token = tokens.shift();
    //     if (token.type()!=KEYWORD || token.text()!="BLUE") ITFAILS;
    // }
  }
  return;
}

//------------------------------------------------------------------------------------------------//

int main(int argc, char *argv[]) {
  rtt_dsxx::ScalarUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    tstConsole_Token_Stream(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstConsole_Token_Stream.cc
//------------------------------------------------------------------------------------------------//
