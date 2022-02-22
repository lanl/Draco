//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   c4/test/tstParallelUnitTest.cc
 * \author Kelly Thompson
 * \date   Thu Jul 7 2011
 * \brief  Test the functionality of the class ParallelUnitTest
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

#include "c4/C4_Functions.hh"
#include "c4/ParallelUnitTest.hh"
#include "ds++/Release.hh"
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
void tstMemberFunctions(ParallelUnitTest &ut, stringstream &output) {
  // test check_all function for failing case. We put this first so we can flush
  // output afterwards.
  {
    string const msg("Testing the check_all member function for failing case.");
    ut.check_all(rtt_c4::node() == 0, msg);

    string const data(output.str());
    size_t const found = data.find(msg);
    if (ut.numPasses == 0 && found != string::npos) {
      cout << "Test: passed\n\t check_all member function works for failing "
              "case ."
           << endl;
      ut.reset();
      output.clear();
    } else {
      cout << "Test: failed\n\t passes member function failed for failing case." << endl;
    }
  }

  // test pass functions
  {
    string const msg("Testing the passes member function.");
    ut.passes(msg);

    string const data(output.str());
    size_t const found = data.find(msg);
    if (ut.numPasses == 1 && found != string::npos)
      cout << "Test: passed\n\t passes member function works." << endl;
    else
      cout << "Test: failed\n\t passes member function failed." << endl;
  }

  // test check_all function
  {
    string const msg("Testing the check_all member function for passing case.");
    ut.check_all(true, msg);

    string const data(output.str());
    size_t const found = data.find(msg);
    if (ut.numPasses == 1 && found != string::npos)
      cout << "Test: passed\n\t check_all member function works for pass." << endl;
    else
      cout << "Test: failed\n\t check_all member function failed." << endl;
  }

  // Test the PASSMSG macro
  {
    string msg("Testing the PASSMSG macro.");
    PASSMSG(msg);

    string const data(output.str());
    size_t const found = data.find(msg);
    if (ut.numPasses == 2 && found != string::npos)
      cout << "Test: passed\n\tPASSMSG macro works." << endl;
    else
      cout << "Test: passed\n\tPASSMSG macro failed." << endl;
  }

  // Try the status member function
  {
    ostringstream msg;
    msg << "Done testing tstParallelUnitTest on node " << node() << ".";
    string const expectedString(msg.str());
    ut.status();

    string const data(output.str());
    size_t const found = data.find(expectedString);
    if (ut.numPasses == 2 && found != string::npos)
      cout << "Test: passed\n\t status() member function works." << endl;
    else
      cout << "Test: failed\n\t status() member function failed." << endl;
  }

  return;
}

//------------------------------------------------------------------------------------------------//

int main(int argc, char *argv[]) {
  int retVal(0);
  stringstream output; // capture output here.

  try {
    {
      rtt_c4::ParallelUnitTest ut(argc, argv, release, output);
      tstMemberFunctions(ut, output);
    } // closing scope should call the destructor for ParallelUnitTest

    // Since we are capturing the output in a stringstream, we must echo the
    // output to stdout so that ctest can pick up the 'passed' message.
    cout << output.str() << endl;
  } catch (rtt_dsxx::assertion &err) {
    if (err.what() == string("Success")) // expected value for --verion
                                         // cmd line option
    {
      // Since we are capturing the output in a stringstream, we must echo the
      // output to stdout so that ctest can pick up the 'passed' message.
      cout << output.str() << endl;
    } else {
      cout << "ERROR: While testing " << argv[0] << ", " << err.what() << endl;
      retVal++;
    }
  } catch (exception &err) {
    cout << "ERROR: While testing " << argv[0] << ", " << err.what() << endl;
    retVal++;
  } catch (...) {
    cout << "ERROR: While testing " << argv[0] << ", "
         << "An unknown exception was thrown." << endl;
    retVal++;
  }

  return retVal;
}

//------------------------------------------------------------------------------------------------//
// end of tstunit_test.cc
//------------------------------------------------------------------------------------------------//
