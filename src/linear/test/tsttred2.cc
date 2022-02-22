//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   linear/test/tsttred2.cc
 * \author Kent Budge
 * \date   Mon Aug  9 13:06:56 2004
 * \brief  Test the tred2 nonlinear equation solver.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"
#include "ds++/Soft_Equivalence.hh"
#include "linear/tred2.hh"

using namespace std;
using namespace rtt_dsxx;
using namespace rtt_linear;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//
void tsttred2(UnitTest &ut) {
  vector<double> A(4 * 4);
  for (unsigned i = 0; i < 4; i++) {
    for (unsigned j = 0; j < 4; j++) {
      A[i + 4 * j] = i + j + 1;
      if (i == j)
        A[i + 4 * j] += 2;
    }
  }
  vector<double> A0 = A;

  vector<double> d, e;
  tred2(A, 4, d, e);

  vector<double> H(4 * 4, 0.0);
  for (unsigned i = 0; i < 4; i++) {
    H[i + 4 * i] = d[i];
    if (i > 0)
      H[i - 1 + 4 * i] = e[i];
    if (i < 3)
      H[i + 1 + 4 * i] = e[i + 1];
  }

  vector<double> QH(4 * 4);
  for (unsigned i = 0; i < 4; i++) {
    for (unsigned j = 0; j < 4; ++j) {
      double sum = 0.0;
      for (unsigned k = 0; k < 4; ++k) {
        sum += A[i + 4 * k] * H[k + 4 * j];
      }
      QH[i + 4 * j] = sum;
    }
  }
  unsigned err = 0;
  for (unsigned i = 0; i < 4; i++) {
    for (unsigned j = 0; j < 4; ++j) {
      double sum = 0.0;
      for (unsigned k = 0; k < 4; ++k) {
        sum += QH[i + 4 * k] * A[j + 4 * k];
      }
      if (!soft_equiv(sum, A0[i + 4 * j])) {
        err++;
      }
    }
  }
  if (err) {
    ut.failure("tred2 is NOT correct");
  } else {
    ut.passes("tred2 is correct");
  }
}

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  ScalarUnitTest ut(argc, argv, release);
  try {
    tsttred2(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tsttred2.cc
//------------------------------------------------------------------------------------------------//
