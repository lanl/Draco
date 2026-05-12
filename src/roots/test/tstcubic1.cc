//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   roots/test/tstcubic1.cc
 * \author Kent G. Budge
 * \note   Copyright (C) 2010-2026 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "dsxx/Release.hh"
#include "dsxx/ScalarUnitTest.hh"
#include "dsxx/Soft_Equivalence.hh"
#include "dsxx/UnitTest.hh"
#include "roots/cubic1.hh"

namespace {

using namespace std;
using namespace rtt_dsxx;
using namespace rtt_roots;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

void tstcubic1(UnitTest &ut) {
  // Solve (x*x+1)*(x-1) = x*x*x-x*x+x-1 = 0

  double const root = cubic1(-1., 1., -1.);

  if (soft_equiv(root, 1.0)) {
    ut.passes("Correctly solved cubic equation");
  } else {
    ut.failure("Did NOT correctly solve cubic equation");
  }
}

} // namespace

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  rtt_dsxx::ScalarUnitTest ut(argc, argv, release);
  try {
    tstcubic1(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstcubic1.cc
//------------------------------------------------------------------------------------------------//
