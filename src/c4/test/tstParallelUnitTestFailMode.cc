//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   c4/test/tstParallelUnitTestFailMode.cc
 * \author Kelly Thompson
 * \date   Thu Jun 1 17:42:58 2006
 * \brief  Test the functionality of the class ParallelUnitTest
 * \note   Copyright (C) 2011-2022 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

/*
 * NOTE: This is not a good example of how ParallelUnitTest should be used.
 * For an example see c4/tstParallelUnitTest.cc
 *
 * This unit test is setup to check parts of ParallelUnitTest that should not
 * be exposed during normal use.
 */

#include "c4/C4_Functions.hh"
#include "c4/ParallelUnitTest.hh"
#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"
#include "ds++/SystemCall.hh"
#include <sstream>

using namespace std;
using namespace rtt_dsxx;
using namespace rtt_c4;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

/*
 * Demonstrate that the normal access functions work as indended.
 */
void tstOne(UnitTest &ut, std::stringstream &output) {
  bool pass(true);

  // test fail functions:
  {
    string msg("Testing the failure member function.");
    ut.failure(msg);

    string const data(output.str());
    size_t found = data.find(msg);
    if (ut.numFails == 1 && found != string::npos) {
      pass &= true;
      cout << "Test: passed\n\t failure member function works." << endl;
    } else {
      pass &= false;
      cout << "Test: failed\n\t failure member function failed." << endl;
    }
  }

  {
    string msg("Testing the FAILMSG macro.");
    FAILMSG(msg);

    string const data(output.str());
    size_t found = data.find(msg);
    if (ut.numFails == 2 && found != string::npos) {
      pass &= true;
      cout << "Test: passed\n\tFAILMSG macro works." << endl;
    } else {
      pass &= false;
      cout << "Test: passed\n\tFAILMSG macro failed." << endl;
    }
  }

  {
    string msg("Test: failed on line");
    ITFAILS;

    string const data(output.str());
    size_t found = data.find(msg);
    if (ut.numFails == 3 && found != string::npos) {
      pass &= true;
      cout << "Test: passed\n\tITFAILS macro works." << endl;
    } else {
      pass &= false;
      cout << "Test: passed\n\tITFAILS macro failed." << endl;
    }
  }

  // If all is good so far clear the fail messages/counts
  if (pass) {
    cout << "Clearing data from ParalleUnitTest..." << endl;
    ut.reset();
    output.str(""); // empty the data
    output.clear(); // reset the flags
  }

  // test pass functions
  {
    string msg("Testing the passes member function.");
    ut.passes(msg);

    string const data(output.str());
    size_t found = data.find(msg);
    if (ut.numPasses == 1 && found != string::npos)
      cout << "Test: passed\n\t passes member function works." << endl;
    else
      cout << "Test: failed\n\t passes member function failed." << endl;
  }

  {
    string msg("Testing the PASSMSG macro.");
    PASSMSG(msg);

    string const data(output.str());
    size_t found = data.find(msg);
    if (ut.numPasses == 2 && found != string::npos)
      cout << "Test: passed\n\tPASSMSG macro works." << endl;
    else
      cout << "Test: passed\n\tPASSMSG macro failed." << endl;
  }

  // keep the data in the output since a failure really is a failure.

  return;
}

void tstTwo(UnitTest &ut) {
  global_barrier();

  // Only issue a failure on PE 1
  if (node() == 1)
    ut.failure("Forced failure on PE 1");
  else
    ut.passes("Pass on PE != 1");

  return;
}

//------------------------------------------------------------------------------------------------//

int main(int argc, char *argv[]) {
  stringstream output;
  bool node_is_0(false);
  bool numnodes_gt_1(false);
  bool tstOnePass(false);
  bool tstTwoPass(false);
  try {
    // Test ctor for ParallelUnitTest (also tests UnitTest ctor and member
    // function setTestName).
    ParallelUnitTest ut(argc, argv, release, output);
    tstOne(ut, output);
    // save the status
    tstOnePass = (ut.numPasses > 0) && (ut.numFails == 0);

    // Test the case where a failure is only recorded on 1 PE.
    node_is_0 = (node() == 0);
    if (nodes() > 1) {
      tstTwo(ut);
      numnodes_gt_1 = true;
    } else {
      tstTwoPass = true; // assume pass for 1 PE job.
    }

    // ut will destruct now.  Part of this process is a global sum of
    // numPasses and numFails.  We must check these values after the try
    // block.
  } catch (rtt_dsxx::assertion &err) {
    std::string msg = err.what();
    if (msg != std::string("Success")) {
      cout << "ERROR: DBC assertion while testing " << argv[0] << ", " << err.what() << endl;
      // ut.numFails++;
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

  //------------------------------------------------------------------------------------------------//
  // Now that ut has called ~ParallelUnitTest(), we need to examine the
  // output to show that the global sum ran correctly.

  string msg(output.str());
  int retval(0);
  if (numnodes_gt_1 && node_is_0) {
    string const expectedOutput =
        string("tstParallelUnitTestFailMode") + rtt_dsxx::exeExtension + string(" Test: FAILED.");
    size_t found = msg.find(expectedOutput);
    if (found != string::npos)
      tstTwoPass = true;
    else
      tstTwoPass = false;
  }

  // Overall test status:
  if (node_is_0) {
    if (tstOnePass && tstTwoPass) {
      cout << "*********************************************\n"
           << "**** tstParallelUnitTestFailMode Test: PASSED.\n"
           << "*********************************************\n"
           << endl;
      retval = 0;
    } else {
      cout << "*********************************************\n"
           << "**** tstParallelUnitTestFailMode Test: FAILED.\n"
           << "*********************************************\n"
           << endl;
      retval = 1;
    }
  }

  return retval;
}

//------------------------------------------------------------------------------------------------//
// end of tstunit_test.cc
//------------------------------------------------------------------------------------------------//
