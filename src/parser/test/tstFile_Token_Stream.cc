//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   parser/test/tstFile_Token_Stream.cc
 * \author Kent G. Budge
 * \date   Feb 18 2003
 * \brief  Unit tests for File_Token_Stream class.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

#include "c4/ParallelUnitTest.hh"
#include "ds++/Release.hh"
#include "parser/File_Token_Stream.hh"
#include <sstream>

using namespace std;
using namespace rtt_parser;
using namespace rtt_dsxx;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

void tstFile_Token_Stream(rtt_dsxx::UnitTest &ut) {
  // Build path for the input file "scanner_test.inp"
  string const inputFile(ut.getTestSourcePath() + std::string("scanner_test.inp"));

  {
    File_Token_Stream tokens(inputFile);
    tokens.comment("begin tests of File_Token_Stream");
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
    tokens.check_semantics(false, "dummy error");
    tokens.check_semantics(true, "dummy error");
    if (tokens.error_count() != 2)
      FAILMSG("Dummy error NOT counted properly");
    else
      PASSMSG("Dummy error counted properly");

    try {
      throw invalid_argument("dummy exception");
    } catch (exception &msg) {
      tokens.report_semantic_error(msg);
    }
    if (tokens.error_count() != 3)
      FAILMSG("Dummy exception NOT reported properly");
    else
      PASSMSG("Dummy exception reported properly");

    tokens.open(inputFile);

    token = tokens.lookahead(4);
    if (token.type() != KEYWORD || token.text() != "BLACK")
      FAILMSG("lookahead(4) does NOT have correct value after open");
    else
      PASSMSG("lookahead(4) has correct value  after open");
  }

  {
    set<char> ws;
    ws.insert(':');
    File_Token_Stream tokens(inputFile, ws);
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

    try {
      tokens.check_syntax(true, "dummy syntax error");
      PASSMSG("Syntax error correctly checked");
    } catch (const Syntax_Error & /*msg*/) {
      FAILMSG("Syntax error NOT correctly checked");
    }
    try {
      tokens.check_syntax(false, "dummy syntax error");
      FAILMSG("Syntax error NOT correctly checked");
    } catch (const Syntax_Error & /*msg*/) {
      PASSMSG("Syntax error correctly checked");
    }
    if (tokens.error_count() != 2) {
      FAILMSG("Syntax errors NOT correctly counted");
    } else {
      PASSMSG("Syntax errors correctly counted");
    }

    token = tokens.shift();
    if (token.type() != KEYWORD || token.text() != "COLOR")
      ITFAILS;

    token = tokens.shift();
    if (token.type() != OTHER || token.text() != "=")
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

    token = tokens.shift();
    if (token.type() != EXIT)
      ITFAILS;
    token = tokens.shift();
    if (token.type() != EXIT)
      ITFAILS;

    tokens.rewind();
    token = tokens.lookahead();
    token = tokens.shift();
    if (token.type() != KEYWORD || token.text() != "BLUE")
      ITFAILS;
  }

  //------------------------------------------------------------------------------------------------//

  {
    try {
      File_Token_Stream tokens("no such file");
      // The preceeding file does not exist.
      ostringstream errmsg;
      errmsg << "File_Token_Stream did not throw an expected exception.\n"
             << "\tThe constructor should throw an exception if the requested\n"
             << "\tfile can not be opened." << endl;
      FAILMSG(errmsg.str());
    } catch (invalid_argument const & /*a*/) {
      std::ostringstream errmsg;
      errmsg << "File_Token_Stream threw an expected exception.\n"
             << "\tThe constructor should throw an exception if the requested\n"
             << "\tfile can not be opened." << endl;
      PASSMSG(errmsg.str());
    } catch (...) {
      ostringstream errmsg;
      errmsg << "File_Token_Stream threw an unknown exception "
             << "during contruction." << endl;
      FAILMSG(errmsg.str());
    }
  }
  {
    try {
      set<char> ws;
      File_Token_Stream tokens("no such file", ws);
      // The preceeding file does not exist.
      ostringstream errmsg;
      errmsg << "File_Token_Stream did not throw an expected exception.\n"
             << "\tThe constructor should throw an exception if the requested\n"
             << "\tfile can not be opened." << endl;
      FAILMSG(errmsg.str());
    } catch (invalid_argument const & /*a*/) {
      std::ostringstream errmsg;
      errmsg << "File_Token_Stream threw an expected exception.\n"
             << "\tThe constructor should throw an exception if the requested\n"
             << "\tfile can not be opened." << endl;
      PASSMSG(errmsg.str());
    } catch (...) {
      ostringstream errmsg;
      errmsg << "File_Token_Stream threw an unknown exception "
             << "during contruction." << endl;
      FAILMSG(errmsg.str());
    }
  }
  {
    try {
      File_Token_Stream tokens;
      tokens.open("no such file");
      // The preceeding file does not exist.
      ostringstream errmsg;
      errmsg << "File_Token_Stream did not throw an expected exception.\n"
             << "\tThe constructor should throw an exception if the requested\n"
             << "\tfile can not be opened." << endl;
      FAILMSG(errmsg.str());
    } catch (invalid_argument const & /*a*/) {
      std::ostringstream errmsg;
      errmsg << "File_Token_Stream threw an expected exception.\n"
             << "\tThe constructor should throw an exception if the requested\n"
             << "\tfile can not be opened." << endl;
      PASSMSG(errmsg.str());
    } catch (...) {
      ostringstream errmsg;
      errmsg << "File_Token_Stream threw an unknown exception "
             << "during contruction." << endl;
      FAILMSG(errmsg.str());
    }
  }

  //-------------------------------------------------------------------------//
  {
    // Build path for the input file "scanner_recovery.inp"
    string const inputFile2(ut.getTestSourcePath() + std::string("scanner_recovery.inp"));

    File_Token_Stream tokens;
    tokens.open(inputFile2);
    // bool exception = false;
    try {
      tokens.shift();
      ostringstream msg;
      msg << "Token_Stream did not throw an exception when\n"
          << "\tunbalanced quotes were read from the input\n"
          << "\tfile, \"scanner_recover.inp\" (line 1)." << endl;
      FAILMSG(msg.str());
    } catch (const Syntax_Error &msg) {
      string errmsg = msg.what();
      string expected("syntax error");
      if (errmsg == expected) {
        ostringstream message;
        message << "Caught expected exception from Token_Stream.\n"
                << "\tunbalanced quotes were read from the input\n"
                << "\tfile, \"scanner_recover.inp\" (line 1)." << endl;
        PASSMSG(message.str());
      } else
        ITFAILS;
    }

    try {
      tokens.shift();
      ostringstream msg;
      msg << "Token_Stream did not throw an exception when\n"
          << "\tunbalanced quotes were read from the input\n"
          << "\tfile, \"scanner_recover.inp\" (line 2)." << endl;
      FAILMSG(msg.str());
    } catch (const Syntax_Error &msg) {
      //cout << msg.what() << endl;
      // exception = true;
      string errmsg = msg.what();
      string expected("syntax error");
      if (errmsg == expected) {
        ostringstream message;
        message << "Caught expected exception from Token_Stream.\n"
                << "\tunbalanced quotes were read from the input\n"
                << "\tfile, \"scanner_recover.inp\" (line 2)." << endl;
        PASSMSG(message.str());
      } else
        ITFAILS;
    }
  }

  // Test #include directive.
  {
    File_Token_Stream tokens(ut.getTestSourcePath() + std::string("parallel_include_test.inp"));

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
    File_Token_Stream dummy;
    ut.check(dummy.lookahead().type() == EXIT, "empty stream returns EXIT");
  }
  return;
}

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  rtt_c4::ParallelUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    Insist(rtt_c4::nodes() == 1, "This test requires exactly 1 PE.");
    tstFile_Token_Stream(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstFile_Token_Stream.cc
//------------------------------------------------------------------------------------------------//
