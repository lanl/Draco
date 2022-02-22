//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   parser/test/tstParallel_File_Token_Stream.cc
 * \author Kent Budge
 * \date   Fri Apr  4 09:34:28 2003
 * \brief  Unit tests for class Parallel_File_Token_Stream
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

#include "c4/ParallelUnitTest.hh"
#include "ds++/Release.hh"
#include "parser/Parallel_File_Token_Stream.hh"
#include <cmath>
#include <sstream>

using namespace std;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

void tstParallel_File_Token_Stream(rtt_dsxx::UnitTest &ut) {
  using namespace rtt_parser;

  // Build path for the input file "scanner_test.inp"
  string const inputFile(ut.getTestSourcePath() + std::string("scanner_test.inp"));

  {
    Parallel_File_Token_Stream tokens(inputFile);
    tokens.comment("begin test of Parallel_File_Token_Stream");
    if (tokens.whitespace() != Text_Token_Stream::default_whitespace)
      FAILMSG("Whitespace not set correctly");
    else
      PASSMSG("Whitespace set correctly.");

    Token token = tokens.lookahead(4);
    if (token.type() != KEYWORD || token.text() != "BLACK")
      FAILMSG("Keyword not read correctly");
    else
      PASSMSG("Keyword read correctly.");

    tokens.report_semantic_error(token, "dummy error");
    if (tokens.error_count() != 1)
      FAILMSG("Semantic error not handled correctly.");
    else
      PASSMSG("Semantic error handled correctly.");

    tokens.report_semantic_error("dummy error");
    if (tokens.error_count() != 2)
      FAILMSG("Second semantic error not handled correctly.");
    else
      PASSMSG("Second semantic error handled correctly.");
  }

  {
    set<char> ws;
    ws.insert(':');
    Parallel_File_Token_Stream tokens(inputFile, ws);
    if (tokens.whitespace() != ws)
      FAILMSG("Whitespace not set correctly");
    else
      PASSMSG("Whitespace set correctly.");

    Token token = tokens.lookahead(4);
    if (token.type() != OTHER || token.text() != "=")
      FAILMSG("'=' token not read correctly");
    else
      PASSMSG("'=' token read correctly.");

    token = tokens.shift();
    if (token.type() != KEYWORD || token.text() != "BLUE")
      FAILMSG("Keyword BLUE not read correctly");
    else
      PASSMSG("Keyword BLUE read correctly.");

    token = tokens.lookahead();
    if (token.type() != KEYWORD || token.text() != "GENERATE ERROR")
      FAILMSG("Keyword GENERATE ERROR not read correctly");
    else
      PASSMSG("Keyword GENERATE ERROR read correctly.");

    token = tokens.shift();
    if (token.type() != KEYWORD || token.text() != "GENERATE ERROR")
      ITFAILS;

    token = tokens.shift();
    if (token.type() != KEYWORD || token.text() != "GENERATE ANOTHER ERROR")
      ITFAILS;

    token = Token('$', "test_parser");
    tokens.pushback(token);

    token = tokens.shift();
    if (token.type() != OTHER || token.text() != "$")
      ITFAILS;

    bool caught(false);
    try {
      tokens.report_syntax_error(token, "dummy syntax error");
    } catch (const Syntax_Error &) {
      caught = true;
      {
        ostringstream msg;
        msg << "Parallel_File_Token_Stream threw an expected exception when\n"
            << "\ta syntax error was reported by Token_Stream." << endl;
        PASSMSG(msg.str());
      }
    }
    FAIL_IF_NOT(caught);

    FAIL_IF(tokens.error_count() != 1);

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

    // Check invariance even when --with-dbc=0.
    FAIL_IF(!tokens.check_class_invariants());
  }

  {
    try {
      Parallel_File_Token_Stream tokens("no such file");
      // The preceeding file does not exist.
      ostringstream errmsg;
      errmsg << "Parallel_File_Token_Stream did not throw an expected exception.\n"
             << "\tThe constructor should throw an exception if the requested\n"
             << "\tfile can not be opened." << endl;
      FAILMSG(errmsg.str());
      // Token token = tokens.shift();
      // if (token.type()!=ERROR) ITFAILS;
    } catch (std::invalid_argument const & /*a*/) {
      std::ostringstream errmsg;
      errmsg << "Parallel_File_Token_Stream threw an expected exception.\n"
             << "\tThe constructor should throw an exception if the requested\n"
             << "\tfile can not be opened." << endl;
      PASSMSG(errmsg.str());
    } catch (...) {
      ostringstream errmsg;
      errmsg << "Parallel_File_Token_Stream threw an unknown exception "
             << "during contruction." << endl;
      FAILMSG(errmsg.str());
    }
  }

  {
    // Build path for the input file "scanner_test.inp"
    string const inputFile2(ut.getTestSourcePath() + std::string("scanner_recovery.inp"));

    Parallel_File_Token_Stream tokens(inputFile2);
    bool exception = false;
    try {
      tokens.shift();
    } catch (const Syntax_Error &msg) {
      cout << msg.what() << endl;
      exception = true;
    }
    if (!exception)
      ITFAILS;

    exception = false;
    try {
      tokens.shift();
    } catch (const Syntax_Error &msg) {
      cout << msg.what() << endl;
      exception = true;
    }
    if (!exception)
      ITFAILS;

    // Try reopening
    string const inputFile3(ut.getTestSourcePath() + std::string("scanner_test.inp"));

    tokens.open(inputFile3);

    Token token = tokens.lookahead(4);
    if (token.type() != KEYWORD || token.text() != "BLACK") {
      FAILMSG("Keyword not read correctly");
    } else {
      PASSMSG("Keyword read correctly.");
    }
  }

  {
    // Test default constructor
    string const inputFile3(ut.getTestSourcePath() + std::string("scanner_test.inp"));

    Parallel_File_Token_Stream tokens;
    tokens.open(inputFile3);

    Token token = tokens.lookahead(4);
    if (token.type() != KEYWORD || token.text() != "BLACK") {
      FAILMSG("Keyword not read correctly");
    } else {
      PASSMSG("Keyword read correctly.");
    }
  }

  // Test #include directive.
  {
    Parallel_File_Token_Stream tokens(ut.getTestSourcePath() +
                                      std::string("parallel_include_test.inp"));

    Token token = tokens.shift();
    ut.check(token.text() == "topmost", "parse top file in include sequence");
    token = tokens.shift();
    ut.check(token.text() == "second", "parse included file in include sequence");
    token = tokens.shift();
    ut.check(token.text() == "topmost2", "parse top file after include sequence");

    // Try rewind
    tokens.rewind();
    token = tokens.shift();
    ut.check(token.text() == "topmost", "parse top file in include sequence");
    token = tokens.shift();
    ut.check(token.text() == "second", "parse included file in include sequence");
    token = tokens.shift();
    ut.check(token.text() == "topmost2", "parse top file after include sequence");

    // Try open of file in middle of include
    tokens.rewind();
    token = tokens.shift();
    ut.check(token.text() == "topmost", "parse top file in include sequence");
    token = tokens.shift();
    ut.check(token.text() == "second", "parse included file in include sequence");
    tokens.open(ut.getTestSourcePath() + std::string("parallel_include_test.inp"));
    token = tokens.shift();
    ut.check(token.text() == "topmost", "parse top file in include sequence");

    // Try rewind in middle of include
    tokens.rewind();
    token = tokens.shift();
    ut.check(token.text() == "topmost", "parse top file in include sequence");
    token = tokens.shift();
    ut.check(token.text() == "second", "parse included file in include sequence");
    tokens.rewind();
    token = tokens.shift();
    ut.check(token.text() == "topmost", "parse top file in include sequence");

    // Check empty stream
    Parallel_File_Token_Stream dummy;
    ut.check(dummy.lookahead().type() == EXIT, "empty stream returns EXIT");
  }
}

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  rtt_c4::ParallelUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    tstParallel_File_Token_Stream(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstParallel_File_Token_Stream.cc
//------------------------------------------------------------------------------------------------//
