//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   linear/test/tstqrdcmp.cc
 * \author Kent Budge
 * \date   Mon Aug  9 13:39:20 2004
 * \brief  Unit teests for qrdcmp
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"
#include "ds++/Soft_Equivalence.hh"
#include "linear/qrdcmp.hh"
#include <array>
#include <iostream>
#include <vector>

using namespace std;
using namespace rtt_dsxx;
using namespace rtt_linear;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//
void tstqrdcmp(UnitTest &ut) {
  vector<double> A = {2., 1., 3., 5.};
  vector<double> C, D;
  qrdcmp(A, 2, C, D);

  // Compute QR to verify
  array<double, 2> uj = {A[0 + 2 * 0], A[1 + 2 * 0]};
  array<double, 4> Qj = {1 - uj[0] * uj[0] / C[0], -uj[1] * uj[0] / C[0], -uj[0] * uj[1] / C[0],
                         1 - uj[1] * uj[1] / C[0]};

  array<double, 4> QR = {Qj[0 + 2 * 0] * D[0], Qj[1 + 2 * 0] * D[0],
                         Qj[0 + 2 * 0] * A[0 + 2 * 1] + Qj[0 + 2 * 1] * D[1],
                         Qj[1 + 2 * 0] * A[0 + 2 * 1] + Qj[1 + 2 * 1] * D[1]};

  if (soft_equiv(QR[0 + 2 * 0], 2.0)) {
    PASSMSG("0,0 is correct");
  } else {
    FAILMSG("0,0 is NOT correct");
  }
  if (soft_equiv(QR[0 + 2 * 1], 3.0)) {
    PASSMSG("0,1 is correct");
  } else {
    FAILMSG("0,1 is NOT correct");
  }
  if (soft_equiv(QR[1 + 2 * 0], 1.0)) {
    PASSMSG("1,0 is correct");
  } else {
    FAILMSG("1,0 is NOT correct");
  }
  if (soft_equiv(QR[1 + 2 * 1], 5.0)) {
    PASSMSG("1,1 is correct");
  } else {
    FAILMSG("1,1 is NOT correct");
  }
}

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  ScalarUnitTest ut(argc, argv, release);
  try {
    tstqrdcmp(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstqrdcmp.cc
//------------------------------------------------------------------------------------------------//
