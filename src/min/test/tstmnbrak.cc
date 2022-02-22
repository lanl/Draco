//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   min/test/tstmnbrak.cc
 * \author Kent Budge
 * \date   Tue Aug 26 13:12:30 2008
 * \brief  
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

//------------------------------------------------------------------------------------------------//

#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"
#include "min/mnbrak.hh"

using namespace std;
using namespace rtt_dsxx;
using namespace rtt_min;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

double func(double const x) { return x * x; }

double func3(double const x) { return (x + 1) * x * (x - 1); }

//------------------------------------------------------------------------------------------------//
void tstmnbrak(UnitTest &ut) {
  double ax, bx, cx, fa, fb, fc;

  ax = 1.0;
  bx = 2.0;

  mnbrak(ax, bx, cx, fa, fb, fc, func);

  if (min(ax, min(bx, cx)) < 0.0 && max(ax, max(bx, cx)) > 0.0) {
    ut.passes("minimum bracketed");
  } else {
    ut.failure("minimum NOT bracketed");
  }
}

//------------------------------------------------------------------------------------------------//

int main(int argc, char *argv[]) {
  ScalarUnitTest ut(argc, argv, release);
  try {
    tstmnbrak(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstmnbrak.cc
//------------------------------------------------------------------------------------------------//
