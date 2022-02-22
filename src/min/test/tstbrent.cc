//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   min/test/tstbrent.cc
 * \author Kent G. Budge
 * \date   Tue Nov 16 17:26:03 2010
 * \brief  
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

//------------------------------------------------------------------------------------------------//

#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"
#include "ds++/Soft_Equivalence.hh"
#include "min/brent.hh"

using namespace std;
using namespace rtt_dsxx;
using namespace rtt_min;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

double f(double x) { return cos(x); }

void tstbrent(UnitTest &ut) {
  double xmin;
  double root = brent(0.0, 6.28, 2.0, f, 1.0e-12, xmin);

  cout << xmin << " " << root << endl;

  if (soft_equiv(xmin, M_PI)) {
    ut.passes("correctly found first minimum of cos");
  } else {
    ut.failure("did NOT correctly find first minimum of cos");
  }
}

//------------------------------------------------------------------------------------------------//

int main(int argc, char *argv[]) {
  ScalarUnitTest ut(argc, argv, release);
  try {
    tstbrent(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstbrent.cc
//------------------------------------------------------------------------------------------------//
