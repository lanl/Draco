//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   ds++/test/tstCheck_Strings.cc
 * \author John M. McGhee
 * \date   Sun Jan 30 14:57:09 2000
 * \brief  Test code for the Check_Strings utility functions.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "ds++/Check_Strings.hh"
#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"

using namespace std;

//------------------------------------------------------------------------------------------------//
void Check_Strings_test(rtt_dsxx::UnitTest &ut) {

  // Define a vector of strings for testing
  vector<string> names = {"this", "is",  "a#", "test", "xxx!", "space check", "123", "x",
                          "test", "dog", "is", "cat",  "",     "abc"};

  // Print a header
  cout << "\n*** String Utilities Test Program ***\n\n";

  // List the test string

  cout << "The " << names.size() << " strings to be tested are: \n";
  for (const auto &name : names)
    cout << "\"" << name << "\"\n";
  cout << endl;

  //----------------------------------------------------------------------------------------------//
  // Test for illegal characters.

  cout << "Illegal character utility test:" << endl;
  string const bad_chars = "()[]* !^#$/";
  auto result = rtt_dsxx::check_string_chars(names.begin(), names.end(), bad_chars);
  if (result.size() == 0) {
    FAILMSG("Failed to find bad characters in string definition.");
  } else {
    PASSMSG("Successfully found bad characters in string definition.");
    for (const auto &bad_entry : result)
      cout << "Found disallowed character(s) in string: \"" << *bad_entry << "\"" << endl;
    cout << "The following characters are forbidden:\n"
         << " \"" << bad_chars << "\","
         << " as well as any white-space characters." << endl;
  }

  if (result.size() == 3)
    PASSMSG("result.size() == 3");
  if (ut.numFails == 0) {
    FAIL_IF_NOT(*result[0] == "a#" && *result[1] == "xxx!" && *result[2] == "space check");
  }

  if (ut.numFails == 0)
    PASSMSG("*** Illegal character function test: PASSED ***");
  else
    FAILMSG("*** Illegal character function test: FAILED ***");

  //----------------------------------------------------------------------------------------------//
  // Test for acceptable lengths.

  cout << "String length utility test:" << endl;
  int const low = 1;
  int const high = 4;
  auto result2 = rtt_dsxx::check_string_lengths(names.begin(), names.end(), low, high);
  if (result2.size() == 0) {
    FAILMSG("Failed to find bad characters in string definition.");
  } else {
    PASSMSG("Successfully found bad characters in string definition.");
    for (const auto &bad_entry : result2)
      cout << "Size of string is not in allowable range: \"" << *bad_entry << "\"" << endl;
    cout << "Strings lengths must be greater than " << low << " and less than " << high << "."
         << endl;
  }

  FAIL_IF_NOT(result2.size() == 2);
  if (ut.numFails == 0) {
    FAIL_IF_NOT(*result2[0] == "space check" && *result2[1] == "");
  }
  if (ut.numFails == 0)
    PASSMSG("*** String length function test: PASSED ***");
  else
    FAILMSG("*** String length function test: FAILED ***");

  //----------------------------------------------------------------------------------------------//
  // Test for unique names.

  cout << "Unique strings utility test:" << endl;
  auto result3 = rtt_dsxx::check_strings_unique(names.begin(), names.end());
  if (result3.size() == 0) {
    FAILMSG("Failed to find bad characters in string definition.");
  } else {
    PASSMSG("Successfully found bad characters in string definition.");
    for (const auto &bad_entry : result3)
      cout << "Duplicate string found: \"" << *bad_entry << "\"\n";
    cout << "All strings must be unique!" << endl;
  }

  FAIL_IF_NOT(result3.size() == 2);
  if (ut.numFails == 0) {
    FAIL_IF_NOT(*result3[0] == "is" && *result3[1] == "test");
  }
  if (ut.numFails == 0)
    PASSMSG("*** Unique string function test: PASSED ***");
  else
    FAILMSG("*** Unique string function test: FAILED ***");

  return;
}

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  rtt_dsxx::ScalarUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    Check_Strings_test(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstCheck_Strings.cc
//------------------------------------------------------------------------------------------------//
