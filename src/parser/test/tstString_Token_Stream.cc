//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   parser/test/tstString_Token_Stream.cc
 * \author Kent G. Budge
 * \date   Feb 18 2003
 * \brief  Unit tests for String_Token_Stream class.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"
#include "parser/String_Token_Stream.hh"
#include "parser/utilities.hh"
#include <sstream>

using namespace std;
using namespace rtt_parser;
using namespace rtt_dsxx;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

void tstString_Token_Stream(UnitTest &ut) {
  // Build path for the input file "scanner_test.inp"
  string const stInputFile(ut.getTestSourcePath() + std::string("scanner_test.inp"));

  ifstream infile(stInputFile.c_str());
  string contents;
  while (true) {
    char const c = static_cast<char>(infile.get());
    if (infile.eof() || infile.fail())
      break;
    contents += c;
  }

  {
    String_Token_Stream tokens(contents);
    tokens.comment("begin test of String_Token_Stream");
    if (tokens.whitespace() != Text_Token_Stream::default_whitespace)
      FAILMSG("whitespace characters are NOT correct defaults");
    else
      PASSMSG("whitespace characters are correct defaults");

    Token token = tokens.lookahead(4);
    if (token.type() != KEYWORD || token.text() != "BLACK")
      FAILMSG("lookahead(4) does NOT have correct value");
    else
      PASSMSG("lookahead(4) has correct value");

    tokens.report_semantic_error(token, "dummy error");
    if (tokens.error_count() != 1)
      FAILMSG("Dummy error NOT counted properly");
    else
      PASSMSG("Dummy error counted properly");

    if (!tokens.check_class_invariants())
      ITFAILS;
  }

  {
    set<char> ws({':'});
    String_Token_Stream tokens(contents, ws);
    if (tokens.whitespace() != ws)
      FAILMSG("whitespace characters are NOT correctly specified");
    else
      PASSMSG("whitespace characters are correctly specified");

    Token token = tokens.lookahead(4);
    if (token.type() != OTHER || token.text() != "=")
      FAILMSG("lookahead(4) does NOT have correct value");
    else
      PASSMSG("lookahead(4) has correct value");

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

    if (tokens.error_count() != 1) {
      FAILMSG("Syntax error NOT correctly counted");
    } else {
      PASSMSG("Syntax error correctly counted");
      if (tokens.messages() == "test_parser\ndummy syntax error\n")
        PASSMSG("Correct error message");
      else
        FAILMSG("NOT correct error message");
    }

    token = tokens.shift();
    FAIL_IF(token.type() != KEYWORD || token.text() != "COLOR");

    token = tokens.shift();
    FAIL_IF(token.type() != OTHER || token.text() != "=");

    token = tokens.shift();
    FAIL_IF(token.type() != KEYWORD || token.text() != "BLACK");

    token = tokens.shift();
    FAIL_IF(token.type() != END);

    token = tokens.shift();
    FAIL_IF(token.type() != OTHER || token.text() != "-");

    token = tokens.shift();
    FAIL_IF(token.type() != REAL || token.text() != "1.563e+3");

    token = tokens.shift();
    FAIL_IF(token.type() != REAL || token.text() != "1.563e+3");

    token = tokens.shift();
    FAIL_IF(token.type() != REAL || token.text() != ".563e+3");

    token = tokens.shift();
    FAIL_IF(token.type() != OTHER || token.text() != ".");

    token = tokens.shift();
    FAIL_IF(token.type() != OTHER || token.text() != "-");

    token = tokens.shift();
    FAIL_IF(token.type() != REAL || token.text() != "1.");

    token = tokens.shift();
    FAIL_IF(token.type() != REAL || token.text() != "1.563");

    token = tokens.shift();
    FAIL_IF(token.type() != REAL || token.text() != "1.e+3");

    token = tokens.shift();
    FAIL_IF(token.type() != REAL || token.text() != "1.e3");

    token = tokens.shift();
    FAIL_IF(token.type() != REAL || token.text() != "1e+3");

    token = tokens.shift();
    FAIL_IF(token.type() != REAL || token.text() != "1e3");

    token = tokens.shift();
    FAIL_IF(token.type() != INTEGER || token.text() != "19090");

    token = tokens.shift();
    FAIL_IF(token.type() != INTEGER || token.text() != "01723");

    token = tokens.shift();
    FAIL_IF(token.type() != INTEGER || token.text() != "0x1111a");

    token = tokens.shift();
    FAIL_IF(token.type() != INTEGER || token.text() != "0");

    token = tokens.shift();
    FAIL_IF(token.type() != INTEGER || token.text() != "8123");

    token = tokens.shift();
    FAIL_IF(token.type() != STRING || token.text() != "\"manifest string\"");

    token = tokens.shift();
    FAIL_IF(token.type() != STRING || token.text() != R"("manifest \"string\"")");

    token = tokens.shift();
    FAIL_IF(token.type() != OTHER || token.text() != "@");

    token = tokens.shift();
    FAIL_IF(token.type() != INTEGER || token.text() != "1");

    token = tokens.shift();
    FAIL_IF(token.type() != KEYWORD || token.text() != "e");

    token = tokens.shift();
    FAIL_IF(token.type() != INTEGER || token.text() != "0");

    token = tokens.shift();
    FAIL_IF(token.type() != KEYWORD || token.text() != "x");

    token = tokens.shift();
    FAIL_IF(token.type() != EXIT);
    token = tokens.shift();
    FAIL_IF(token.type() != EXIT);

    tokens.rewind();
    token = tokens.lookahead();
    token = tokens.shift();
    FAIL_IF(token.type() != KEYWORD || token.text() != "BLUE");
  }

  //-------------------------------------------------------------------------//
  {
    // Build path for the input file "scanner_recovery.inp"
    string const srInputFile(ut.getTestSourcePath() + std::string("scanner_recovery.inp"));

    ifstream linfile(srInputFile.c_str());
    string lcontents;
    while (true) {
      char const c = static_cast<char>(linfile.get());
      if (linfile.eof() || linfile.fail())
        break;
      lcontents += c;
    }
    String_Token_Stream tokens(lcontents);
    try {
      tokens.shift();
      ostringstream msg;
      msg << "Token_Stream did not throw an exception when unbalanced quotes were read from the \n"
          << "\tinput file, \"scanner_recover.inp\" (line 1)." << endl;
      FAILMSG(msg.str());
    } catch (const Syntax_Error &msg) {
      // cout << msg.what() << endl;
      // exception = true;
      string errmsg = msg.what();
      string expected("syntax error");
      if (errmsg == expected) {
        ostringstream message;
        message << "Caught expected exception from Token_Stream. unbalanced quotes were read \n"
                << "\tfrom the input file, \"scanner_recover.inp\" (line 1)." << endl;
        PASSMSG(message.str());
      } else
        ITFAILS;
    }

    try {
      tokens.shift();
      ostringstream msg;
      msg << "Token_Stream did not throw an exception when unbalanced quotes were read from the \n"
          << "\tinput file, \"scanner_recover.inp\" (line 2)." << endl;
      FAILMSG(msg.str());
    } catch (const Syntax_Error &msg) {
      //cout << msg.what() << endl;
      // exception = true;
      string errmsg = msg.what();
      string expected("syntax error");
      if (errmsg == expected) {
        ostringstream message;
        message << "Caught expected exception from Token_Stream. unbalanced quotes were read \n"
                << "\tfrom the input file, \"scanner_recover.inp\" (line 2)." << endl;
        PASSMSG(message.str());
      } else
        ITFAILS;
    }

    // Now test assignment of new string and dipthong OTHER tokens

    String_Token_Stream tokens2 = string("<= >= && ||");

    Token token = tokens2.shift();
    FAIL_IF(token.text() != "<=");

    token = tokens2.shift();
    FAIL_IF(token.text() != ">=");

    token = tokens2.shift();
    FAIL_IF(token.text() != "&&");

    token = tokens2.shift();
    FAIL_IF(token.text() != "||");

    token = tokens2.shift();
    FAIL_IF(token.type() != EXIT);
  }

  {
    String_Token_Stream tokens("09");
    if (tokens.is_nb_whitespace('\t')) {
      PASSMSG("tab correctly identified as nonbreaking whitespace");
    } else {
      FAILMSG("tab NOT correctly identified as nonbreaking whitespace");
    }
    Token token = tokens.shift();
    if (token.type() != INTEGER || token.text() != "0") {
      FAILMSG("did NOT scan 09 correctly");
    } else {
      PASSMSG("scanned 09 correctly");
    }
  }
  {
    String_Token_Stream tokens("_, __, _ _, > < & | 1E3 0XA");
    FAIL_IF(tokens.shift().text() != "_");
    FAIL_IF(tokens.shift().text() != "__");
    FAIL_IF(tokens.shift().text() != "_ _");
    FAIL_IF(tokens.shift().text() != ">");
    FAIL_IF(tokens.shift().text() != "<");
    FAIL_IF(tokens.shift().text() != "&");
    FAIL_IF(tokens.shift().text() != "|");
    FAIL_IF(!soft_equiv(parse_real(tokens), 1e3));
    FAIL_IF(parse_integer(tokens) != 10);
  }

  // Test that missing closing quote is a syntax error.
  {
    String_Token_Stream tokens("\"quote");
    try {
      tokens.shift();
      FAILMSG("Did NOT correctly report missing closing quote as syntax error");
    } catch (const Syntax_Error & /*msg*/) {
      PASSMSG("missing closing quote correctly thrown and caught");
    }
  }

  // Test that #include is treated as a syntax error (not supported)
  {
    String_Token_Stream tokens("#include \"dummy.inp\"");
    try {
      tokens.shift();
      FAILMSG("Did NOT correctly report #include as error");
    } catch (const Syntax_Error & /*msg*/) {
      cout << "expected: " << tokens.messages() << endl;
      PASSMSG("#include not supported error correctly thrown and caught");
    }
  }

  // Test that # without a valid directive is treated as syntax error
  {
    String_Token_Stream tokens("# !");
    try {
      tokens.shift();
      FAILMSG("Did NOT correctly report #! as error");
    } catch (const Syntax_Error & /*msg*/) {
      PASSMSG("invalid #directive correctly thrown and caught");
    }
  }
  {
    String_Token_Stream tokens("#bad");
    try {
      tokens.shift();
      FAILMSG("Did NOT correctly report #bad as error");
    } catch (const Syntax_Error & /*msg*/) {
      PASSMSG("invalid #bad correctly thrown and caught");
    }
  }
  {
    String_Token_Stream tokens("#include, bad");
    try {
      tokens.shift();
      FAILMSG("Did NOT correctly report #include, bad as error");
    } catch (const Syntax_Error & /*msg*/) {
      PASSMSG("invalid #insist, bad correctly thrown and caught");
    }
  }

  return;
}

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  ScalarUnitTest ut(argc, argv, release);
  try {
    tstString_Token_Stream(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstString_Token_Stream.cc
//------------------------------------------------------------------------------------------------//
