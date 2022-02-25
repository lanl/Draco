//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   min/test/tstpowell.cc
 * \author Kent Budge
 * \date   Tue Jul  6 10:00:38 2004
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"
#include "min/powell.hh"
#include "units/PhysicalConstants.hh"

using namespace std;
using namespace rtt_min;
using namespace rtt_dsxx;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

unsigned constexpr NP = 2;

double xf(double const x, vector<double> const &a) {
  double const xh = sqrt(x);
  double Result = (sqrt(8.0) + a[0] * xh + a[1] * x + x * xh) * x * xh;
  Result = fabs(Result / cube(sqrt(x * x + 2 * x)) - 1);
  return Result;
}

double func(vector<double> const &a) {
  double Result = 0.0;
  double x(1.0e-5);
  while (x < 1000.0) {
    Result = max(Result, xf(x, a));
    x *= 1.1;
  }
  return Result;
}

void tstpowell(UnitTest &ut) {
  vector<double> p(NP, 0.0);
  vector<double> xi(NP * NP, 0.0);
  for (unsigned i = 0; i < NP; ++i) {
    xi[i + NP * i] = 1.0;
  }
  unsigned iter = 10000;
  double fret(0);
  double tolerance(1.0e-5);

  powell(p, xi, tolerance, iter, fret, func);

  for (unsigned i = 0; i < NP; ++i)
    cout << "a[" << i << "] = " << p[i] << endl;

  cout << "Maximum error: " << fret << endl;

  std::array<double, 2> tmp = {1.34601, 4.19265e-09};
  vector<double> expectedSolution(tmp.begin(), tmp.end());
  if (rtt_dsxx::soft_equiv(p.begin(), p.end(), expectedSolution.begin(), expectedSolution.end(),
                           tolerance)) {
    ut.passes("Found expected solution.");
  } else {
    ut.failure("Did not find expected solution.");
  }
  return;
}

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  ScalarUnitTest ut(argc, argv, release);
  try {
    tstpowell(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of testpowell.cc
//------------------------------------------------------------------------------------------------//
