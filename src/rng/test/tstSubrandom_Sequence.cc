//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   rng/test/tstSubrandom_Sequence.cc
 * \author Kent Budge
 * \date   Thu Dec 22 14:16:45 2006
 * \brief  Test the Subrandom_Sequence class
 * \note   Copyright (C) 2006-2022 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"
#include "ds++/Soft_Equivalence.hh"
#include "rng/LC_Subrandom_Generator.hh"
#include <fstream>
#include <iostream>

using namespace std;
using namespace rtt_dsxx;
using namespace rtt_rng;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

void tstSubrandom_Sequence(UnitTest &ut) {
  {
    // Test LC_Subrandom_Generator

    LC_Subrandom_Generator lcsg;
    double value = lcsg.shift();
    FAIL_IF_NOT(soft_equiv(value, 0.999741748906672));
    lcsg.shift_vector();
  }

  if (ut.numFails == 0)
    PASSMSG("test subrandom passed");
  return;
}

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  ScalarUnitTest ut(argc, argv, release);
  try {
    tstSubrandom_Sequence(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstSubrandom_Sequence.cc
//------------------------------------------------------------------------------------------------//
