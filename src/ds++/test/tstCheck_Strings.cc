//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ds++/test/tstCheck_Strings.cc
 * \author John M. McGhee
 * \date   Sun Jan 30 14:57:09 2000
 * \brief  Test code for the Check_Strings utility functions.
 * \note   Copyright (C) 2016 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "ds++/Check_Strings.hh"
#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"

using namespace std;

//---------------------------------------------------------------------------//

void Check_Strings_test(rtt_dsxx::UnitTest &ut) {
  // Define a vector of strings for testing
  string n[] = {"this", "is",  "a#", "test", "xxx!", "space check", "123", "x",
                "test", "dog", "is", "cat",  "",     "abc"};
  const int nn = sizeof(n) / sizeof(string);
  vector<string> names(&n[0], &n[nn]);
  typedef vector<string>::iterator VS_iter;

  // Print a header
  cout << "\n*** String Utilities Test Program ***\n" << endl;

  // List the test string

  cout << "The " << names.size() << " strings to be tested are: " << endl;
  for (size_t i = 0; i < names.size(); ++i)
    cout << "\"" << names[i] << "\"" << endl;
  cout << endl;

  //---------------------------------------------------------------------------//

  // Test for illegal characters.

  cout << "Illegal character utility test:" << endl;
  string bad_chars = "()[]* !^#$/";
  vector<VS_iter> result =
      rtt_dsxx::check_string_chars(names.begin(), names.end(), bad_chars);
  if (result.size() == 0) {
    FAILMSG("Failed to find bad characters in string definition.");
  } else {
    PASSMSG("Successfully found bad characters in string definition.");
    for (size_t i = 0; i < result.size(); i++)
      cout << "Found disallowed character(s) in string: \"" << *result[i]
           << "\"" << endl;
    cout << "The following characters are forbidden:" << endl
         << " \"" << bad_chars << "\","
         << " as well as any white-space characters." << endl;
  }

  if (result.size() == 3)
    PASSMSG("result.size() == 3");
  if (ut.numFails == 0) {
    if (*result[0] == "a#" && *result[1] == "xxx!" &&
        *result[2] == "space check")
      PASSMSG("result has expected values.");
    else
      FAILMSG("result did not have expected values.");
  }

  if (ut.numFails == 0)
    PASSMSG("*** Illegal character function test: PASSED ***");
  else
    FAILMSG("*** Illegal character function test: FAILED ***");

  //---------------------------------------------------------------------------//

  // Test for acceptable lengths.

  cout << "String length utility test:" << endl;
  int low = 1;
  int high = 4;
  vector<VS_iter> result2 =
      rtt_dsxx::check_string_lengths(names.begin(), names.end(), low, high);
  if (result2.size() == 0) {
    FAILMSG("Failed to find bad characters in string definition.");
  } else {
    PASSMSG("Successfully found bad characters in string definition.");
    for (size_t i = 0; i < result2.size(); i++)
      cout << "Size of string is not in allowable range: \"" << *result2[i]
           << "\"" << endl;
    cout << "Strings lengths must be greater than " << low << " and less than "
         << high << "." << endl;
  }

  if (result2.size() == 2)
    PASSMSG("result2 has the expected size.");
  if (ut.numFails == 0) {
    if (*result2[0] == "space check" && *result2[1] == "")
      PASSMSG(" result2 has expected content.");
    else
      FAILMSG(" result2 does not have the expected content.");
  }
  if (ut.numFails == 0)
    PASSMSG("*** String length function test: PASSED ***");
  else
    FAILMSG("*** String length function test: FAILED ***");

  //---------------------------------------------------------------------------//

  // Test for unique names.

  cout << "Unique strings utility test:" << endl;
  vector<VS_iter> result3 =
      rtt_dsxx::check_strings_unique(names.begin(), names.end());
  if (result3.size() == 0) {
    FAILMSG("Failed to find bad characters in string definition.");
  } else {
    PASSMSG("Successfully found bad characters in string definition.");
    for (size_t i = 0; i < result3.size(); i++)
      cout << "Duplicate string found: \"" << *result3[i] << "\"" << endl;
    cout << "All strings must be unique!" << endl;
  }

  if (result3.size() != 2)
    ITFAILS;
  if (ut.numFails == 0) {
    if (*result3[0] == "is" && *result3[1] == "test")
      PASSMSG("result3 has expected content.");
    else
      FAILMSG("result3 did not have the expected content.");
  }
  if (ut.numFails == 0)
    PASSMSG("*** Unique string function test: PASSED ***");
  else
    FAILMSG("*** Unique string function test: FAILED ***");

  return;
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[]) {
  rtt_dsxx::ScalarUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    Check_Strings_test(ut);
  }
  UT_EPILOG(ut);
}

//---------------------------------------------------------------------------//
// end of tstCheck_Strings.cc
//---------------------------------------------------------------------------//
