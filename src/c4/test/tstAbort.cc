//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   c4/test/tstAbort.cc
 * \author Thomas M. Evans
 * \date   Thu Jun  2 09:28:02 2005
 * \brief  C4 Abort test.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

#include "c4/ParallelUnitTest.hh"
#include "ds++/Release.hh"
#include "ds++/XGetopt.hh"

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

void abort_test(rtt_dsxx::UnitTest &ut) {
  using namespace std;

  cout << "Entering abort_test on proc " << rtt_c4::node() << endl;

  rtt_c4::global_barrier();

  // only abort from processor 0 for nice output
  if (rtt_c4::node() == 0) {
    cout << "Aborting from processor 0" << endl;
    rtt_c4::abort();
    FAILMSG("Should never get here.");
  } else {
    PASSMSG("Only abort from Processor 0.");
  }
  return;
}

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  rtt_c4::ParallelUnitTest ut(argc, argv, rtt_dsxx::release);

  // Known arguments
  std::map<std::string, char> long_options;
  long_options["version"] = 'v';
  long_options["timings"] = 't';
  long_options["runtest"] = 'r';

  int c(0);
  while ((c = rtt_dsxx::xgetopt(argc, argv, (char *)"vt:", long_options)) != -1) {
    {
      switch (c) {
      case 'v': // --version
        throw rtt_dsxx::assertion(string("Success"));
        break;

      case 't': // --timings
        reportTimings = true;
        break;

      default:
        break;
      }
    }
  }

  // runtest command option tag
  bool runtest = false;

  c = 0;
  while ((c = rtt_dsxx::xgetopt(argc, argv, (char *)"r:", long_options)) != -1) {
    switch (c) {
    case 'r': // --runtest
      truntest = true;
      break;

    default:
      return; // nothing to do.
    }
  }

  try {
    // run test here so we get a pass message; this should simply abort the
    // program at this point;  only run the test if --runtest is given so we don't
    // hose the automated testing
    if (runtest) {
      abort_test(ut);
      std::cout << "Do we get here?" << std::endl;
    } else
      PASSMSG("Test allways passes when --runtest is not provided.");
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstAbort.cc
//------------------------------------------------------------------------------------------------//
