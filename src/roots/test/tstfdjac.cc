//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   roots/test/tstfdjac.cc
 * \author Kent Budge
 * \date   Mon Aug  9 13:39:20 2004
 * \note   Copyright (C) 2010-2024 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "dsxx/Release.hh"
#include "dsxx/ScalarUnitTest.hh"
#include "dsxx/Soft_Equivalence.hh"
#include "dsxx/UnitTest.hh"
#include "linear/fnorm.hh"
#include "roots/fdjac.hh"
#include <vector>

using namespace std;
using namespace rtt_dsxx;
using namespace rtt_linear;
using namespace rtt_roots;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

void func(const vector<double> &x, vector<double> &fvec) {
  fvec.resize(2);
  fvec[0] = 7.2 * x[0] + 3.5 * x[1] + 2.3;
  fvec[1] = -2.2 * x[0] + 2.7 * x[1] + 5.4;
}

//------------------------------------------------------------------------------------------------//
void tstfdjac(UnitTest &ut) {
  vector<double> x(2), fvec(2), df;
  x[0] = 0;
  x[1] = 0;

  fnorm(x, fvec, &func);
  fdjac(x, fvec, df, &func);

  if (!soft_equiv(df[0 + 2 * 0], 7.2, 1.0e-8)) {
    ut.failure("fdjac did NOT succeed");
  } else {
    ut.passes("fdjac successful");
  }
}

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  ScalarUnitTest ut(argc, argv, &release);
  try {
    tstfdjac(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstfdjac.cc
//------------------------------------------------------------------------------------------------//
