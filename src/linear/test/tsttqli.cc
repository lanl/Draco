//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   linear/test/tsttqli.cc
 * \author Kent Budge
 * \date   Mon Aug  9 13:06:56 2004
 * \brief  Test the tqli eigenvector solver.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"
#include "ds++/Soft_Equivalence.hh"
#include "linear/tqli.hh"
#include "linear/tred2.hh"

using namespace std;
using namespace rtt_dsxx;
using namespace rtt_linear;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

//------------------------------------------------------------------------------------------------//
void tsttqli(UnitTest &ut) {
  vector<double> A(2 * 2);

  A[0 + 2 * 0] = 1.8;
  A[0 + 2 * 1] = 2.2;
  A[1 + 2 * 0] = 2.2;
  A[1 + 2 * 1] = 3.2;

  vector<double> d, e;

  // Trivial example
  tred2(A, 2, d, e);
  tqli(d, e, 2, A);

  double e0 = d[0];
  double xx0 = A[0 + 2 * 0];
  double xx1 = A[1 + 2 * 0];

  double x0 = (1.8 - e0) * xx0 + 2.2 * xx1;
  double x1 = 2.2 * xx0 + (3.2 - e0) * xx1;

  if (fabs(x0) < 1.0e-12 && fabs(x1) < 1.0e-12) {
    ut.passes("first eigenvector is correct");
  } else {
    ut.failure("first eigenvector is NOT correct");
  }

  double e1 = d[1];
  xx0 = A[0 + 2 * 1];
  xx1 = A[1 + 2 * 1];

  x0 = (1.8 - e1) * xx0 + 2.2 * xx1;
  x1 = 2.2 * xx0 + (3.2 - e1) * xx1;

  if (fabs(x0) < 1.0e-12 && fabs(x1) < 1.0e-12) {
    ut.passes("second eigenvector is correct");
  } else {
    ut.failure("second eigenvector is NOT correct");
  }
}

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  ScalarUnitTest ut(argc, argv, release);
  try {
    tsttqli(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tsttqli.cc
//------------------------------------------------------------------------------------------------//
