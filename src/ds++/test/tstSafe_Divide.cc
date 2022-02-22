//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   ds++/test/tstSafe_Divide.cc
 * \author Mike Buksas
 * \date   Tue Jun 21 16:02:52 2005
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

#include "ds++/Release.hh"
#include "ds++/Safe_Divide.hh"
#include "ds++/ScalarUnitTest.hh"
#include "ds++/Soft_Equivalence.hh"

using namespace std;
using namespace rtt_dsxx;

#if defined(MSVC)
#pragma warning(push)
// warning C4756: overflow in constant arithmetic
#pragma warning(disable : 4756)
#endif

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//
void test(rtt_dsxx::UnitTest &ut) {
  double const max = numeric_limits<double>::max();

  double const big = 1.0e200;
  double const tiny = 1.0e-200;

  if (!rtt_dsxx::soft_equiv(safe_pos_divide(big, tiny), max))
    ITFAILS;
  if (!rtt_dsxx::soft_equiv(safe_pos_divide(10.0, 5.0), 2.0))
    ITFAILS;

  if (!rtt_dsxx::soft_equiv(safe_divide(big, tiny), max))
    ITFAILS;
  if (!rtt_dsxx::soft_equiv(safe_divide(-big, tiny), -max))
    ITFAILS;
  if (!rtt_dsxx::soft_equiv(safe_divide(-big, -tiny), max))
    ITFAILS;
  if (!rtt_dsxx::soft_equiv(safe_divide(big, -tiny), -max))
    ITFAILS;

  if (!rtt_dsxx::soft_equiv(safe_divide(10.0, 5.0), 2.0))
    ITFAILS;
  if (!rtt_dsxx::soft_equiv(safe_divide(-10.0, 5.0), -2.0))
    ITFAILS;
  if (!rtt_dsxx::soft_equiv(safe_divide(-10.0, -5.0), 2.0))
    ITFAILS;
  if (!rtt_dsxx::soft_equiv(safe_divide(10.0, -5.0), -2.0))
    ITFAILS;

  if (ut.numFails == 0)
    PASSMSG("done with test().");
  return;
}

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  rtt_dsxx::ScalarUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    test(ut);
  }
  UT_EPILOG(ut)
}

#if defined(MSVC)
#pragma warning(pop)
#endif

//------------------------------------------------------------------------------------------------//
//  end of tstSafe_Divide.cc
//------------------------------------------------------------------------------------------------//
