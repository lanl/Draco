//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   linear/test/tstqr_unpack.cc
 * \author Kent Budge
 * \date   Mon Aug  9 13:39:20 2004
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"
#include "ds++/Soft_Equivalence.hh"
#include "linear/qr_unpack.hh"
#include "linear/qrdcmp.hh"

using namespace std;
using namespace rtt_dsxx;
using namespace rtt_linear;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//
void tstqr_unpack(UnitTest &ut) {
  vector<double> A(2 * 2);
  A[0 + 2 * 0] = 2.;
  A[0 + 2 * 1] = 3.;
  A[1 + 2 * 0] = 1.;
  A[1 + 2 * 1] = 5.;

  vector<double> C, D;

  qrdcmp(A, 2, C, D);

  vector<double> Qt;
  vector<double> R = A;
  qr_unpack(R, 2, C, D, Qt);

  if (soft_equiv(R[0 + 2 * 0], D[0]))
    PASSMSG("R[0+2*0] is correct");
  else
    FAILMSG("R[0+2*0] is NOT correct");
}

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  ScalarUnitTest ut(argc, argv, release);
  try {
    tstqr_unpack(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstqr_unpack.cc
//------------------------------------------------------------------------------------------------//
