//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   rng/test/tstSubrandom_Sequence.cc
 * \author Kent Budge
 * \brief  Test the Subrandom_Sequence class
 * \note   Copyright (C) 2006-2026 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "dsxx/Release.hh"
#include "dsxx/ScalarUnitTest.hh"
#include "dsxx/Soft_Equivalence.hh"
#include "dsxx/UnitTest.hh"
#include "rng/LC_Subrandom_Generator.hh"

// cSpell:ignore lcsg

namespace {

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
    double const value = lcsg.shift();
    FAIL_IF_NOT(soft_equiv(value, 0.999741748906672));
    lcsg.shift_vector();
  }

  if (ut.numFails == 0)
    PASSMSG("test subrandom passed");
  return;
}

} // namespace

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
