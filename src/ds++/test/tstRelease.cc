//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ds++/test/tstRelease.cc
 * \author Kelly Thompson <kgt@lanl.gov>
 * \date   Friday, Jul 29, 2016, 10:05 am
 * \brief  Check basic functionality of Release.hh/cc files.
 * \note   Copyright (C) 2011-2016 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//

#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"
#include "ds++/Soft_Equivalence.hh"
#include <iostream>
#include <map>
#include <sstream>
#include <string>

using namespace std;
using namespace rtt_dsxx;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

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
    std::map<std::string, unsigned> wc =
        ut.get_word_count(releaseString, verbose);

    if (wc[string("DRACO_DIAGNOSTICS")] != 1)
      ITFAILS;
    if (wc[string("build")] != 2)
      ITFAILS;
  }

  {
    // ----------------------------------------
    // Print the copyright statement and author list
    ostringstream const copyrightString(copyright());
    cout << "\ncopyright() = \n" << copyrightString.str() << endl;

    if (copyrightString.str().length() > 0)
      PASSMSG("copyrightString len > 0");
    else
      FAILMSG("copyrightString len == 0");

    bool verbose(false);
    std::map<std::string, unsigned> wc =
        ut.get_word_count(copyrightString, verbose);

    if (wc[string("CCS-2")] != 1)
      ITFAILS;
    if (wc[string("Copyright")] != 1)
      ITFAILS;
    if (wc[string("Contributers")] != 1)
      ITFAILS;
    if (wc[string("Team")] != 1)
      ITFAILS;
  }

  return;
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[]) {
  ScalarUnitTest ut(argc, argv, release);
  try {
    maintest(ut);
  }
  UT_EPILOG(ut);
}

//---------------------------------------------------------------------------//
// end of tstJayenneRelease.cc
//---------------------------------------------------------------------------//
