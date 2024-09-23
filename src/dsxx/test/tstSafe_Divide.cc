//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   dsxx/test/tstSafe_Divide.cc
 * \author Mike Buksas
 * \date   Tue Jun 21 16:02:52 2005
 * \note   Copyright (C) 2024 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "dsxx/Release.hh"
#include "dsxx/Safe_Divide.hh"
#include "dsxx/ScalarUnitTest.hh"
#include "dsxx/Soft_Equivalence.hh"
#include "dsxx/UnitTest.hh"
#include <limits>

using namespace std;
using namespace rtt_dsxx;

#if defined(MSVC)
#pragma warning(push)
// warning C4756: overflow in constant arithmetic
#pragma warning(disable : 4756)
#endif

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  rtt_dsxx::ScalarUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    double constexpr max = numeric_limits<double>::max();
    double constexpr big = 1.0e200;
    double constexpr tiny = 1.0e-200;

    FAIL_IF_NOT(rtt_dsxx::soft_equiv(safe_pos_divide(big, tiny), max));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(safe_pos_divide(10.0, 5.0), 2.0));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(safe_divide(big, tiny), max));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(safe_divide(-big, tiny), -max));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(safe_divide(-big, -tiny), max));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(safe_divide(big, -tiny), -max));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(safe_divide(10.0, 5.0), 2.0));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(safe_divide(-10.0, 5.0), -2.0));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(safe_divide(-10.0, -5.0), 2.0));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(safe_divide(10.0, -5.0), -2.0));

    if (ut.numFails == 0)
      PASSMSG("done with test().");
  }
  UT_EPILOG(ut)
}

#if defined(MSVC)
#pragma warning(pop)
#endif

//------------------------------------------------------------------------------------------------//
//  end of tstSafe_Divide.cc
//------------------------------------------------------------------------------------------------//
