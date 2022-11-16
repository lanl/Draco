//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   linear/test/tstsvdcmp.cc
 * \author Kent Budge
 * \date   Mon Aug  9 13:39:20 2004
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"
#include "ds++/Soft_Equivalence.hh"
#include "linear/svdcmp.hh"
#include <array>

using namespace std;
using namespace rtt_dsxx;
using namespace rtt_linear;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//
void tstsvdcmp(UnitTest &ut) {
  vector<double> U = {2., 1., 4., 3., 5., 4.};
  vector<double> W, V;

  svdcmp(U, 3, 2, W, V);

  // Compute U*W*Tr(V) to verify

  array<double, 9> WV = {W[0] * V[0 + 2 * 0], W[1] * V[0 + 2 * 1], W[0] * V[1 + 2 * 0],
                         W[1] * V[1 + 2 * 1]};

  array<double, 6> UWV = {U[0 + 3 * 0] * WV[0 + 2 * 0] + U[0 + 3 * 1] * WV[1 + 2 * 0],
                          U[1 + 3 * 0] * WV[0 + 2 * 0] + U[1 + 3 * 1] * WV[1 + 2 * 0],
                          U[2 + 3 * 0] * WV[0 + 2 * 0] + U[2 + 3 * 1] * WV[1 + 2 * 0],
                          U[0 + 3 * 0] * WV[0 + 2 * 1] + U[0 + 3 * 1] * WV[1 + 2 * 1],
                          U[1 + 3 * 0] * WV[0 + 2 * 1] + U[1 + 3 * 1] * WV[1 + 2 * 1],
                          U[2 + 3 * 0] * WV[0 + 2 * 1] + U[2 + 3 * 1] * WV[1 + 2 * 1]};

  if (soft_equiv(UWV[0 + 3 * 0], 2.0))
    PASSMSG("0,0 is correct");
  else
    FAILMSG("0,0 is NOT correct");
  if (soft_equiv(UWV[0 + 3 * 1], 3.0))
    PASSMSG("0,1 is correct");
  else
    FAILMSG("0,1 is NOT correct");
  if (soft_equiv(UWV[1 + 3 * 0], 1.0))
    PASSMSG("1,0 is correct");
  else
    FAILMSG("1,0 is NOT correct");
  if (soft_equiv(UWV[1 + 3 * 1], 5.0))
    PASSMSG("1,1 is correct");
  else
    FAILMSG("1,1 is NOT correct");
  if (soft_equiv(UWV[2 + 3 * 1], 4.0))
    PASSMSG("2,1 is correct");
  else
    FAILMSG("2,1 is NOT correct");
  if (soft_equiv(UWV[2 + 3 * 1], 4.0))
    PASSMSG("2,1 is correct");
  else
    FAILMSG("2,1 is NOT correct");

  // Now decompose transpose, to exercise different code paths

  U[0 + 2 * 0] = 2.;
  U[1 + 2 * 0] = 3.;
  U[0 + 2 * 1] = 1.;
  U[1 + 2 * 1] = 5.;
  U[0 + 2 * 2] = 4.;
  U[1 + 2 * 2] = 4.;

  svdcmp(U, 2, 3, W, V);

  // Compute U*W*Tr(V) to verify

  WV[0 + 3 * 0] = W[0] * V[0 + 3 * 0];
  WV[0 + 3 * 1] = W[0] * V[1 + 3 * 0];
  WV[0 + 3 * 2] = W[0] * V[2 + 3 * 0];
  WV[1 + 3 * 0] = W[1] * V[0 + 3 * 1];
  WV[1 + 3 * 1] = W[1] * V[1 + 3 * 1];
  WV[1 + 3 * 2] = W[1] * V[2 + 3 * 1];
  WV[2 + 3 * 0] = W[2] * V[0 + 3 * 2];
  WV[2 + 3 * 1] = W[2] * V[1 + 3 * 2];
  WV[2 + 3 * 2] = W[2] * V[2 + 3 * 2];

  UWV[0 + 2 * 0] =
      U[0 + 2 * 0] * WV[0 + 3 * 0] + U[0 + 2 * 1] * WV[1 + 3 * 0] + U[0 + 2 * 2] * WV[2 + 3 * 0];
  UWV[0 + 2 * 1] =
      U[0 + 2 * 0] * WV[0 + 3 * 1] + U[0 + 2 * 1] * WV[1 + 3 * 1] + U[0 + 2 * 2] * WV[2 + 3 * 1];
  UWV[0 + 2 * 2] =
      U[0 + 2 * 0] * WV[0 + 3 * 2] + U[0 + 2 * 1] * WV[1 + 3 * 2] + U[0 + 2 * 2] * WV[2 + 3 * 2];
  UWV[1 + 2 * 0] =
      U[1 + 2 * 0] * WV[0 + 3 * 0] + U[1 + 2 * 1] * WV[1 + 3 * 0] + U[1 + 2 * 2] * WV[2 + 3 * 0];
  UWV[1 + 2 * 1] =
      U[1 + 2 * 0] * WV[0 + 3 * 1] + U[1 + 2 * 1] * WV[1 + 3 * 1] + U[1 + 2 * 2] * WV[2 + 3 * 1];
  UWV[1 + 2 * 2] =
      U[1 + 2 * 0] * WV[0 + 3 * 2] + U[1 + 2 * 1] * WV[1 + 3 * 2] + U[1 + 2 * 2] * WV[2 + 3 * 2];

  if (soft_equiv(UWV[0 + 2 * 0], 2.0))
    PASSMSG("0,0 is correct");
  else
    FAILMSG("0,0 is NOT correct");
  if (soft_equiv(UWV[0 + 2 * 1], 1.0))
    PASSMSG("0,1 is correct");
  else
    FAILMSG("0,1 is NOT correct");
  if (soft_equiv(UWV[0 + 2 * 2], 4.0))
    PASSMSG("0,2 is correct");
  else
    FAILMSG("0,2 is NOT correct");
  if (soft_equiv(UWV[1 + 2 * 0], 3.0))
    PASSMSG("1,0 is correct");
  else
    FAILMSG("1,0 is NOT correct");
  if (soft_equiv(UWV[1 + 2 * 1], 5.0))
    PASSMSG("1,1 is correct");
  else
    FAILMSG("1,1 is NOT correct");
  if (soft_equiv(UWV[1 + 2 * 2], 4.0))
    PASSMSG("1,2 is correct");
  else
    FAILMSG("1,2 is NOT correct");
}

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  ScalarUnitTest ut(argc, argv, release);
  try {
    tstsvdcmp(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstsvdcmp.cc
//------------------------------------------------------------------------------------------------//
