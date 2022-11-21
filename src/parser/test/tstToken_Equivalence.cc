//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   parser/test/tstToken_Equivalence.cc
 * \author Kelly Thompson
 * \date   Fri Jul 21 09:10:49 2006
 * \brief  Unit test for functions in Token_Equivalence.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"
#include "parser/Token_Equivalence.hh"
#include <sstream>

using namespace std;
using namespace rtt_dsxx;
using namespace rtt_parser;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

void tstOne(UnitTest &ut) {
  // create some data
  ostringstream data;
  data << "METIS decomposition specified\n"
       << "Dump cycle interval defaulting to 1.\n"
       << "Cycle       : 0\n"
       << "Time Step   : 1e-16 s.\n"
       << "Problem Time: 0 s.\n"
       << "error(0): 1     spr: 1\n"
       << "error(1): 0.00272636     spr: 0.00272636\n"
       << "error(2): 8.14886e-06     spr: 0.00298892\n"
       << "pid[0] done error(2): 8.14886e-06  spr: 0.00298892\n"
       << "User Cpu time this time step: 5.17 \n"
       << endl;

  // create a string token stream from the data

  String_Token_Stream tokens(data.str());

  // Test Token_Equivalence functions:

  // look for an int associated with a keyword
  check_token_keyword_value(tokens, "Cycle", 0, ut);

  // look for a double associated with a keyword (4th occurance).
  check_token_keyword_value(tokens, "spr", 0.00298892, ut, 4);

  // look for a keyword.
  check_token_keyword(tokens, "User Cpu time this time step", ut);

  return;
}

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  ScalarUnitTest ut(argc, argv, release);
  try {
    tstOne(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstToken_Equivalence.cc
//------------------------------------------------------------------------------------------------//
