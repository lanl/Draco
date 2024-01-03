//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   special_functions/test/tstPower.cc
 * \author Mike Buksas
 * \date   Mon Jul 24 13:47:58 2006
 * \note   Copyright (C) 2010-2024 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"
#include "ds++/Soft_Equivalence.hh"
#include "special_functions/Power.hh"

using namespace std;
using namespace rtt_dsxx;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//
void test(rtt_dsxx::UnitTest &ut) {
  using rtt_sf::Power;

  // Power of an integer:

  if (Power<12>(2) != 4096)
    FAILMSG("2^12 == 4096 in integers failed");
  else
    PASSMSG("");

  // Floating point bases:

  if (!soft_equiv(Power<4>(2.0), 16.0))
    FAILMSG("2.0^4 = 16.0 in float failed.");

  if (!soft_equiv(Power<17>(1.1), 5.054470284992945))
    FAILMSG("1.1^17 failed.");

  if (Power<0>(1) != 1)
    FAILMSG("1^0 in int failed.");

  if (!soft_equiv(Power<0>(1.0), 1.0))
    FAILMSG("1.0^0 in float failed.");
}

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  rtt_dsxx::ScalarUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    test(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstPower.cc
//------------------------------------------------------------------------------------------------//
