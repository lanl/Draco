//----------------------------------*-c++-*----------------------------------//
/*!
 * \file   ds++/test/tstFile_Streams.cc
 * \author Rob Lowrie
 * \date   Sun Nov 21 19:36:12 2004
 * \brief  Tests File_Input and File_Output.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

#include "ds++/File_Streams.hh"
#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"
#include "ds++/Soft_Equivalence.hh"
#include <sstream>

using namespace std;
using rtt_dsxx::File_Input;
using rtt_dsxx::File_Output;
using rtt_dsxx::soft_equiv;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

void test_fileio(rtt_dsxx::UnitTest &ut, const bool is_fformat_binary) {

  string const filename("file_streams." + string(is_fformat_binary ? "binary" : "ascii"));

  int i = 5;
  string s = "  a string with spaces  ";
  double x = 5.6;
  bool bf = false;
  bool bt = true;

  // write the data

  {
    File_Output f(filename, is_fformat_binary);
    f << i;

    // here's how you write strings:
    size_t ssize = s.size();
    f << ssize;
    for (size_t k = 0; k < ssize; k++)
      f << s[k];

    f << x << bf << bt;
  }

  // read the data and make sure it's the same
  {
    int i_in;
    double x_in;
    string s_in;
    bool bf_in;
    bool bt_in;

    File_Input f(filename);
    f >> i_in;

    FAIL_IF_NOT(i == i_in);

    // here's how you read strings:
    size_t ssize;
    f >> ssize;
    FAIL_IF_NOT(ssize == s.size());
    s_in.resize(ssize);
    for (size_t k = 0; k < ssize; k++)
      f >> s_in[k];

    FAIL_IF_NOT(s == s_in);

    f >> x_in >> bf_in >> bt_in;

    FAIL_IF_NOT(soft_equiv(x, x_in));
    FAIL_IF_NOT(bf == bf_in);
    FAIL_IF_NOT(bt == bt_in);

    File_Input fnull("");
  }

  // test some corner cases
  {
    File_Output f;
    f.close();

    f.open("File_Stream_last_was_char.txt");
    f << 'c';
    f.close();

    File_Input fr("File_Stream_last_was_char.txt");
    char c;
    fr >> c;
    FAIL_IF_NOT(c == 'c');

    fr.open("File_Stream_last_was_char.txt");
    fr >> c;
    FAIL_IF_NOT(c == 'c');

    f.open("File_Stream_last_was_char.txt", false);
    f.open("File_Stream_last_was_char.txt", false);
    f << 'c';
    f.close();
    fr.open("File_Stream_last_was_char.txt");
    fr >> c;
    FAIL_IF_NOT(c == 'c');
  }

  if (ut.numFails == 0) {
    ostringstream m;
    m << "test_fileio(";
    if (is_fformat_binary)
      m << "binary";
    else
      m << "ascii";
    m << ") ok.";
    PASSMSG(m.str());
  }
  return;
}

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  rtt_dsxx::ScalarUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    test_fileio(ut, false); // ascii
    test_fileio(ut, true);  // binary
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstFile_Streams.cc
//------------------------------------------------------------------------------------------------//
