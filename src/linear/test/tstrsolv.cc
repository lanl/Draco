//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   linear/test/tstrsolv.cc
 * \author Kent Budge
 * \note   Copyright (C) 2010-2026 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "dsxx/Release.hh"
#include "dsxx/ScalarUnitTest.hh"
#include "dsxx/Soft_Equivalence.hh"
#include "dsxx/UnitTest.hh"
#include "linear/rsolv.hh"
#include <vector>

namespace {

using namespace std;
using namespace rtt_dsxx;
using namespace rtt_linear;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//
void tstrsolv(UnitTest &ut) {
  vector<double> R(2 * 2), b(2);
  R[0 + 2 * 0] = 5.0;
  R[0 + 2 * 1] = 7.0;
  R[1 + 2 * 0] = 0.0;
  R[1 + 2 * 1] = 3.0;
  b[0] = 1.0;
  b[1] = 6.0;

  rsolv(R, 2, b);

  if (soft_equiv(b[1], 2.0)) {
    PASSMSG("b[1] is correct");
  } else {
    FAILMSG("b[1] is NOT correct");
  }
  if (soft_equiv(5 * b[0] + 7 * b[1], 1.0)) {
    PASSMSG("b[0] is correct");
  } else {
    FAILMSG("b[0] is NOT correct");
  }
}

} // namespace

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  ScalarUnitTest ut(argc, argv, release);
  try {
    tstrsolv(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstrsolv.cc
//------------------------------------------------------------------------------------------------//
