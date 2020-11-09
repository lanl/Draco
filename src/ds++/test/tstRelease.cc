//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   ds++/test/tstRelease.cc
 * \author Kelly Thompson <kgt@lanl.gov>
 * \date   Friday, Jul 29, 2016, 10:05 am
 * \brief  Check basic functionality of Release.hh/cc files.
 * \note   Copyright (C) 2016-2020 Triad National Security, LLC.
 *         All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "ds++/DracoStrings.hh"
#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"
#include "ds++/Soft_Equivalence.hh"
#include <iostream>
#include <map>
#include <sstream>
#include <string>

using namespace std;
using namespace rtt_dsxx;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

void maintest(UnitTest &ut) {
  {
    // ----------------------------------------
    // Print the release information
    ostringstream const releaseString(release());
    cout << "\nrelease() = \n" << releaseString.str() << "\n" << endl;

    if (releaseString.str().length() > 0)
      PASSMSG("releaseString len > 0");
    else
      FAILMSG("releaseString len == 0");

    bool verbose(false);
    std::map<std::string, unsigned> wc = rtt_dsxx::get_word_count(releaseString, verbose);

    FAIL_IF_NOT(wc[string("DRACO_DIAGNOSTICS")] == 1);
    FAIL_IF_NOT(wc[string("build")] == 2);
  }

  {
    // ----------------------------------------
    // Print the copyright statement and author list
    ostringstream const copyrightString(rtt_dsxx::remove_color(copyright()));
    cout << "\ncopyright() = \n" << copyrightString.str() << endl;

    if (copyrightString.str().length() > 0)
      PASSMSG("copyrightString len > 0");
    else
      FAILMSG("copyrightString len == 0");

    bool verbose(false);
    std::map<std::string, unsigned> wc = rtt_dsxx::get_word_count(copyrightString, verbose);

    FAIL_IF_NOT(wc[string("CCS-2")] == 1);
    FAIL_IF_NOT(wc[string("Copyright")] == 1);
    FAIL_IF_NOT(wc[string("Contributors")] == 1);
    FAIL_IF_NOT(wc[string("Team")] == 1);
  }

  return;
}

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  ScalarUnitTest ut(argc, argv, release);
  try {
    maintest(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstRelease.cc
//------------------------------------------------------------------------------------------------//
