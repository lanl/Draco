//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   linear/test/tstludcmp.cc
 * \author Kent Budge
 * \note   Copyright (C) 2010-2024 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "dsxx/Release.hh"
#include "dsxx/ScalarUnitTest.hh"
#include "dsxx/Soft_Equivalence.hh"
#include "dsxx/UnitTest.hh"
#include "linear/ludcmp.hh"
#include <vector>

using namespace std;
using namespace rtt_dsxx;
using namespace rtt_linear;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//
void tstludcmp(UnitTest &ut) {
  vector<double> U(4);
  U[0 + 2 * 0] = 2.;
  U[0 + 2 * 1] = 3.;
  U[1 + 2 * 0] = 1.;
  U[1 + 2 * 1] = 5.;

  vector<unsigned> indx(2);
  double d;

  ludcmp(U, indx, d);

  vector<double> b(2), x;
  b[0] = 1.0;
  b[1] = 2.0;

  lubksb(U, indx, b);

  if (soft_equiv(1.0, b[0] * 2 + b[1] * 3)) {
    ut.passes("0 is correct");
  } else {
    ut.failure("0 is NOT correct");
  }
  if (soft_equiv(2.0, b[0] * 1 + b[1] * 5)) {
    ut.passes("0 is correct");
  } else {
    ut.failure("0 is NOT correct");
  }
}

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  ScalarUnitTest ut(argc, argv, release);
  try {
    tstludcmp(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstsvbksb.cc
//------------------------------------------------------------------------------------------------//
