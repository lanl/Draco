//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   roots/test/tstbroydn.cc
 * \author Kent Budge
 * \date   Mon Aug  9 13:06:56 2004
 * \brief  Test the broydn nonlinear equation solver.
 * \note   Copyright (C) 2010-2024 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "dsxx/Release.hh"
#include "dsxx/ScalarUnitTest.hh"
#include "dsxx/UnitTest.hh"
#include "roots/broydn.hh"
#include <cmath>
#include <vector>

using namespace std;
using namespace rtt_dsxx;
using namespace rtt_roots;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

//------------------------------------------------------------------------------------------------//
//! Functional form
namespace {
void broydn_test_function(const vector<double> &x, vector<double> &fvec) {
  fvec.resize(2);
  fvec[0] = cos(x[0] + x[1]);
  fvec[1] = sin(x[0] - x[1]);
}
} // namespace

//------------------------------------------------------------------------------------------------//
//! Test class for Broydn root finder
class Broydn_Test_Class {
public:
  void operator()(const vector<double> &x, vector<double> &fvec) const {
    broydn_test_function(x, fvec);
  }

  void operator()(const vector<double> &x, vector<double> &fvec, vector<double> &dfvec) const {
    broydn_test_function(x, fvec);
    dfvec.resize(4);
    dfvec[0 + 2 * 0] = -fvec[1];
    dfvec[0 + 2 * 1] = -fvec[1];
    dfvec[1 + 2 * 0] = fvec[0];
    dfvec[1 + 2 * 1] = -fvec[0];
  }
};

//------------------------------------------------------------------------------------------------//
void tstbroydn(UnitTest &ut) {
  vector<double> x = {0.2, -0.1};

  // Check broydn solution method for function
  broydn(x, 100.0, &broydn_test_function, 1.0e-2);

  if (fabs(cos(x[0] + x[1])) > 1.0e-12 || fabs(sin(x[0] - x[1])) > 1.0e-12) {
    FAILMSG("broydn: FAILED");
  } else {
    PASSMSG("broydn: passed");
  }

  // Check broydn solution method for class
  x[0] += 0.2;
  x[1] -= 0.1;
  broydn(x, 100.0, Broydn_Test_Class(), 1.0e-2);

  if (fabs(cos(x[0] + x[1])) > 1.0e-12 || fabs(sin(x[0] - x[1])) > 1.0e-12) {
    FAILMSG("broydn: FAILED");
  } else {
    PASSMSG("broydn: passed");
  }

  // Check broydn solution with analytic derivatives
  x[0] += 0.2;
  x[1] -= 0.1;
  broydn(x, 100.0, Broydn_Test_Class(), Broydn_Test_Class(), 1.0e-2, 0.0);
  if (fabs(cos(x[0] + x[1])) > 1.0e-12 || fabs(sin(x[0] - x[1])) > 1.0e-12) {
    FAILMSG("broydn: FAILED");
  } else {
    PASSMSG("broydn: passed");
  }
}

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  ScalarUnitTest ut(argc, argv, &release);
  try {
    tstbroydn(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstbroydn.cc
//------------------------------------------------------------------------------------------------//
