//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   ds++/test/tstScalarUnitTest.cc
 * \author Kelly Thompson
 * \date   Thu May 18 17:17:24 2006
 * \brief  Unit test for the ds++ classes UnitTest and ScalarUnitTest.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "ds++/DracoStrings.hh"
#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"
#include <algorithm>
#include <cstring>
#include <fstream>
#include <iterator>
#include <sstream>

using namespace std;
using namespace rtt_dsxx;

// Provide old style call to pass/fail macros.  Use object name unitTest for this unit test. (remove
// the defines found in ScalarUnitTest.hh).
#undef PASSMSG
#undef ITFAILS
#undef FAILURE
#undef FAILMSG
#define PASSMSG(a) unitTest.passes(a)
#define ITFAILS unitTest.failure(__LINE__);
#define FAILURE unitTest.failure(__LINE__, __FILE__);
#define FAILMSG(a) unitTest.failure(a);

//------------------------------------------------------------------------------------------------//
// Helper
//------------------------------------------------------------------------------------------------//
char *convert_string_to_char_ptr(std::string const &s) {
  auto *pc = new char[s.length() + 1];
  std::strncpy(pc, s.c_str(), s.length() + 1);
  return pc;
}

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

void tstOne(UnitTest &unitTest) {
  unitTest.passes("Looks like the passes member function is working.");
  unitTest.check(true, "Also for check version.");
  unitTest.check_all(true, "Also for check_all version.");
  PASSMSG("Looks like the PASSMSG macro is working as a member function.");
  UT_CHECK(unitTest, true);
  return;
}

//------------------------------------------------------------------------------------------------//
void tstTwo(UnitTest &unitTest) {
  unitTest.failure("Looks like the failure member function is working.");
  unitTest.check(false, "Also for check version.");
  unitTest.check_all(false, "Also for check_all version.");
  FAILMSG("Looks like the FAILMSG macro is working.");
  ITFAILS;
  FAILURE;
  UT_CHECK(unitTest, false);

  // Kill report of failures
  unitTest.reset();

  // We need at least one pass.
  PASSMSG("Done with tstTwo.");
  return;
}

//------------------------------------------------------------------------------------------------//
void tstTwoCheck(UnitTest &unitTest, ostringstream &msg) {
  bool verbose(true);
  std::ostringstream const msg_nocolor(rtt_dsxx::remove_color(msg.str()));
  map<string, unsigned> word_list(rtt_dsxx::get_word_count(msg_nocolor, verbose));

  // Check the list of occurrences against the expected values
  if (word_list[string("Test")] == 9)
    unitTest.passes("Found 9 occurrences of \"Test\"");
  else
    unitTest.failure("Did not find expected number of occurrences of \"Test\"");

  if (word_list[string("failed")] != 7)
    unitTest.failure("Did not find 7 occurrences of failure.");
  if (word_list[string("FAILMSG")] != 1)
    unitTest.failure("Found 1 occurrence of \"FAILMSG\"");
  if (word_list[string("failure")] != 1)
    unitTest.failure("Found 1 occurrence of \"failure\"");

  if (word_list[string("macro")] == 1)
    unitTest.passes("Found 1 occurrence of \"macro\"");
  else
    unitTest.failure("Did not find expected number of occurrences of \"macro\"");

  if (word_list[string("working")] == 2)
    unitTest.passes("Found 2 occurrences of \"working\"");
  else
    unitTest.failure("Did not find expected number of occurrences of \"working\"");

  return;
}

//------------------------------------------------------------------------------------------------//
void tstGetWordCountFile(UnitTest &unitTest) {
  cout << "\ntstGetWordCountFile...\n" << endl;

  // Generate a text file
  string filename("tstScalarUnitTest.sample.txt");
  ofstream myfile(filename.c_str());
  if (myfile.is_open()) {
    myfile << "This is a text file.\n"
           << "Used by tstScalarUnitTest::tstGetWordCountFile\n\n"
           << "foo bar baz\n"
           << "foo bar\n"
           << "foo\n\n";
    myfile.close();
  }

  // Now read the file and parse the contents:
  map<string, unsigned> word_list(get_word_count(filename, false));

  // Some output
  cout << "The world_list has the following statistics (word, count):\n" << endl;
  for (auto &it : word_list)
    cout << it.first << "\t::\t" << it.second << endl;

  // Spot checks on file contents:
  if (word_list[string("This")] != 1)
    ITFAILS;
  if (word_list[string("foo")] != 3)
    ITFAILS;
  if (word_list[string("bar")] != 2)
    ITFAILS;
  if (word_list[string("baz")] != 1)
    ITFAILS;

  if (unitTest.numFails == 0) {
    cout << endl;
    PASSMSG("Successfully parsed text file and generated word_list");
  }
  return;
}

//------------------------------------------------------------------------------------------------//
void tstdbcsettersandgetters(UnitTest &unitTest, int argc, char **argv) {
  std::cout << "Testing Design-by-Contract setters and getters "
            << "for the UnitTest class..." << std::endl;

  // Silent version.
  ostringstream messages;

  // DBC = 0 (all off)
  {
    ScalarUnitTest foo(argc, argv, release, messages);
    foo.dbcRequire(false);
    foo.dbcCheck(false);
    foo.dbcEnsure(false);
    if (foo.dbcRequire())
      ITFAILS;
    if (foo.dbcCheck())
      ITFAILS;
    if (foo.dbcEnsure())
      ITFAILS;
    if (foo.dbcOn())
      ITFAILS;
  }
  // DBC = 1 (Require only)
  {
    ScalarUnitTest foo(argc, argv, release, messages);
    foo.dbcRequire(true);
    foo.dbcCheck(false);
    foo.dbcEnsure(false);
    if (!foo.dbcRequire())
      ITFAILS;
    if (foo.dbcCheck())
      ITFAILS;
    if (foo.dbcEnsure())
      ITFAILS;
    if (!foo.dbcOn())
      ITFAILS;
  }
  // DBC = 2 (Check only)
  {
    ScalarUnitTest foo(argc, argv, release, messages);
    foo.dbcRequire(false);
    foo.dbcCheck(true);
    foo.dbcEnsure(false);
    if (foo.dbcRequire())
      ITFAILS;
    if (!foo.dbcCheck())
      ITFAILS;
    if (foo.dbcEnsure())
      ITFAILS;
    if (!foo.dbcOn())
      ITFAILS;
  }
  // DBC = 4 (Ensure only)
  {
    ScalarUnitTest foo(argc, argv, release, messages);
    foo.dbcRequire(false);
    foo.dbcCheck(false);
    foo.dbcEnsure(true);
    if (foo.dbcRequire())
      ITFAILS;
    if (foo.dbcCheck())
      ITFAILS;
    if (!foo.dbcEnsure())
      ITFAILS;
    if (!foo.dbcOn())
      ITFAILS;
  }

  if (unitTest.numPasses > 0 && unitTest.numFails == 0)
    PASSMSG("UnitTest Design-by-Contract setters and getters are working.");

  return;
}

//------------------------------------------------------------------------------------------------//
void tstVersion(UnitTest &unitTest, char *test) {
  // Check version construction

  // Initialize the argument list
  std::vector<std::string> vs_arguments = {test, "a", "--version"};
  auto argc = static_cast<int>(vs_arguments.size());

  // Convert to 'char *' so we can then use &vc[0] as type char**
  std::vector<char *> vc;
  std::transform(vs_arguments.begin(), vs_arguments.end(), std::back_inserter(vc),
                 convert_string_to_char_ptr);

  char **argv = &vc[0];
  try {
    std::cout << "Testing ctor... ScalarUnitTest(argc, argv, release)" << std::endl;
    ScalarUnitTest myut(argc, argv, release);
    FAILMSG("version construction NOT correct");
  } catch (assertion &err) {
    if (!strcmp(err.what(), "Success"))
      PASSMSG("version construction correct");
    else
      FAILMSG("version construction NOT correct");
  } catch (...) {
    FAILMSG("version construction NOT correct");
  }

  // clean-up memory
  for (auto &i : vc)
    delete[] i;
  return;
}

//------------------------------------------------------------------------------------------------//
void tstPaths(UnitTest &unitTest, char *test) {

  using std::string;

  std::cout << "Begin testing tstPaths(ut)..." << std::endl;

  // Checkpoint
  size_t const nf = unitTest.numFails;

  // There are 4 member functions of UnitTest that return paths:
  string const testBinaryDir(unitTest.getTestPath());
  string const testName(unitTest.getTestName());
  string const testBinaryInputDir(unitTest.getTestInputPath());
  string const testSourceDir(unitTest.getTestSourcePath());

  // helper data
  string const thisFile(__FILE__);
  string testName_wo_suffix(testName);
  if (testName.substr(testName.length() - 4, 4) == string(".exe"))
    testName_wo_suffix = testName.substr(0, testName.length() - 4);

  // Report current state
  std::cout << "\nThe unitTest system reports the following paths:"
            << "\n\tTest Path         = " << testBinaryDir << "\n\tTest Name         = " << testName
            << "\n\tBinary Input Path = " << testBinaryInputDir
            << "\n\tSource Input Path = " << testSourceDir << "\n"
            << std::endl;

  // Checks
  string const stest = rtt_dsxx::getFilenameComponent(
      rtt_dsxx::getFilenameComponent(string(test), rtt_dsxx::FC_NATIVE), rtt_dsxx::FC_REALPATH);
  FAIL_IF_NOT(stest == testBinaryDir + testName);
  FAIL_IF_NOT(testName == string("tstScalarUnitTest") + rtt_dsxx::exeExtension);
  if (rtt_dsxx::fileExists(thisFile)) {
    // For some generators (ninja), the logic of this test fails because the __FILE__ CPP macro is
    // relative to the top level source directory.
    if (thisFile != testSourceDir + testName_wo_suffix + ".cc") {
      // 2nd chance for case-insensitive file systems
      string const lc_thisFile = string_tolower(thisFile);
      string const lc_gold = string_tolower(testSourceDir + testName_wo_suffix + ".cc");
      FAIL_IF_NOT(lc_thisFile == lc_gold);
    }
  }

  // CMake should provide cmake_install.cmake at testBinaryInputDir.
  FAIL_IF_NOT(rtt_dsxx::fileExists(testBinaryInputDir + "cmake_install.cmake"));

  // If this is a multi-config build tool, examine the value of buildType.
  string buildType = rtt_dsxx::getFilenameComponent(testBinaryDir, rtt_dsxx::FC_NAME);
  if (buildType != string("test")) {
    // trim trailing Windows or Unix slash, if any
    if (buildType.substr(buildType.length() - 1, 1) == string("\\") ||
        buildType.substr(buildType.length() - 1, 1) == string("/"))
      buildType = buildType.substr(0, buildType.length() - 1);
    std::cout << "This appears to be a multi-config build tool like Xcode or "
              << "Visual Studio where build type = " << buildType << "." << std::endl;
    if (buildType != string("Release") && buildType != string("Debug") &&
        buildType != string("DebWithRelInfo") && buildType != string("MinSizeRel"))
      FAILMSG(string("Unexpected build type = ") + buildType);
  }

  if (unitTest.numFails == nf) // no new failures
    PASSMSG("tstPaths completed successfully.");
  else
    FAILMSG("tstPaths did not complete successfully.");

  return;
}

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  try {
    // Test ctor for ScalarUnitTest (also tests UnitTest ctor and member function setTestName).
    ScalarUnitTest ut(argc, argv, release);

    // Try to print the copyright and author list
    std::cout << copyright() << std::endl;

    tstOne(ut);

    // Silent version.
    ostringstream messages;
    ScalarUnitTest sut(argc, argv, release, messages);
    tstTwo(sut);

    tstTwoCheck(ut, messages);
    tstGetWordCountFile(ut);
    tstdbcsettersandgetters(ut, argc, argv);
    tstVersion(ut, argv[0]);

    tstPaths(ut, argv[0]);

    // Silenced version
    ScalarUnitTest ssut(argc, argv, release, messages, false);
    messages.str("");
    ssut.check(true, "this test must pass");
    ut.check(messages.str().size() == 0, "verbose==false is silent");

    // Check the status member function
    ut.status();
  }

  catch (rtt_dsxx::assertion &err) {
    string msg = err.what();
    if (msg != string("Success")) {
      cout << "ERROR: While testing " << argv[0] << ", " << err.what() << endl;
      return 1;
    }
    return 0;
  } catch (exception &err) {
    cout << "ERROR: While testing " << argv[0] << ", " << err.what() << endl;
    return 1;
  } catch (...) {
    cout << "ERROR: While testing " << argv[0] << ", "
         << "An unknown exception was thrown" << endl;
    return 1;
  }

  return 0;
}

//------------------------------------------------------------------------------------------------//
// end of tstScalarUnitTest.cc
//------------------------------------------------------------------------------------------------//
