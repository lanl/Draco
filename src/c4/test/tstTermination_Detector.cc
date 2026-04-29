//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   c4/test/tstTermination_Detector.cc
 * \author Kent Budge
 * \note   Copyright (C) 2010-2026 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "dsxx/Release.hh"
#include "dsxx/UnitTest.hh"
#include "c4/ParallelUnitTest.hh"
#include "c4/Termination_Detector.hh"

namespace {

using namespace std;
using namespace rtt_dsxx;
using namespace rtt_c4;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

void tstTermDet(UnitTest &ut) {
  Termination_Detector td(1);

  td.init();

  td.update_receive_count(0);
  td.update_send_count(1);
  td.update_work_count(2);

  for (unsigned c = 0; c < 5; ++c) {
    if (td.is_terminated())
      FAILMSG("Termination_Detection did NOT detect nontermination.");
    else
      PASSMSG("Termination_Detection detected nontermination.");
  }

  // Will hang if the unit fails.  Unfortunately, there's no other portable way to test.
  td.update_receive_count(1);

  while (!td.is_terminated()) { /* do nothing */
  };

  return;
}

} // namespace

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  ParallelUnitTest ut(argc, argv, release);
  try {
    tstTermDet(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstTermination_Detector.cc
//------------------------------------------------------------------------------------------------//
