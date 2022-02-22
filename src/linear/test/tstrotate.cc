//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   linear/test/tstrotate.cc
 * \author Kent Budge
 * \date   Mon Aug  9 13:39:20 2004
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"
#include "ds++/Soft_Equivalence.hh"
#include "linear/rotate.hh"

using namespace std;
using namespace rtt_dsxx;
using namespace rtt_linear;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//
void tstrotate(UnitTest &ut) {
  vector<double> A(3 * 3, 0.0);
  for (unsigned i = 0; i < 3; i++)
    A[i + 3 * i] = 1.0;
  vector<double> B = A;

  rotate(A, B, 3, 1, 0.5, 0.5);

  if (soft_equiv(A[0 + 3 * 0], 1.0))
    PASSMSG("A[0+3*0] is correct");
  else
    FAILMSG("A[0+3*0] is NOT correct");
  if (soft_equiv(A[1 + 3 * 1], sqrt(0.5)))
    PASSMSG("A[1+3*1] is correct");
  else
    FAILMSG("A[1+3*1] is NOT correct");
}

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  ScalarUnitTest ut(argc, argv, release);
  try {
    tstrotate(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstrotate.cc
//------------------------------------------------------------------------------------------------//
