//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   min/test/tstmnbrak.cc
 * \note   Copyright (C) 2010-2026 Triad National Security, LLC., All rights reserved.
 *
 * The functions used for testing are plotted in the following images:
 * - \image html tstmnbrak-func2.png "Example of func2(x)"
 * - \image html tstmnbrak-func3.png "Example of func3(x)"
 * - \image html tstmnbrak-func4.png "Example of func4(x)"
 */
//------------------------------------------------------------------------------------------------//

#include "dsxx/Release.hh"
#include "dsxx/ScalarUnitTest.hh"
#include "dsxx/Soft_Equivalence.hh"
#include "dsxx/UnitTest.hh"
#include "min/mnbrak.hh"
#include <algorithm>
#include <cmath>
#include <iostream>

// cSpell:ignore mnbrak, tstmnbrak

namespace {

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

double func(double const x) { return x * x; }
double func2(double const x) { return (x - 1.0) * x; }
double func3(double const x) { return (x + 1.0) * x * (x - 1.0); }
double func4(double const x) { return x * x * x - 2.0 * x * x - 2.0 * x - 2.0; }
double func5(double const x) {
  return x * x * x * x - 2.0 * x * x * x - 3.0 * x * x - 3.0 * x - 3.0;
}
double func6(double const x) { return std::sin(x) + 0.1 * x; }

//------------------------------------------------------------------------------------------------//
void tstmnbrak(rtt_dsxx::UnitTest &ut) {
  double ax(1.0);
  double bx(2.0);
  double cx(0.0);
  double fa(0.0);
  double fb(0.0);
  double fc(0.0);

  { // Test the branch "u is between bx and cx" and "Minimum lies between bx and cx"
    std::cout << "Examine func(x) = x^2, for x in [1.0, 2.0]\n";
    rtt_min::mnbrak(ax, bx, cx, fa, fb, fc, func);
    if (std::min({ax, bx, cx}) < 0.0 && std::max({ax, bx, cx}) > 0.0) {
      // minimum should be at x=0,
      PASSMSG("minimum bracketed (func)");
    } else {
      FAILMSG("minimum NOT bracketed (func)");
    }
  }
  { // Test the branch "u is between cx and upper limit"
    std::cout << "\nExamine func2(x) = x*(x-1), for x in [0.0, 0.1]\n";
    ax = 0.0;
    bx = 0.1;
    rtt_min::mnbrak(ax, bx, cx, fa, fb, fc, func2);
    if (fa > fb && fc > fb && ((ax < bx && cx > bx) || (cx < bx && ax < bx))) {
      // minimum should be at x=0.5,
      PASSMSG("minimum bracketed (func2)");
    } else {
      FAILMSG("minimum NOT bracketed (func2)");
    }
  }
  { // Test the branch "Minimum lies between ax and u"
    std::cout << "\nExamine func4(x) = x^3 - 2x^2 - 2x - 2, for x in [1.5, 3.5]\n";
    ax = 1.5;
    bx = 3.5;
    rtt_min::mnbrak(ax, bx, cx, fa, fb, fc, func4);
    if ((fa > fb && fc > fb) && ((ax > bx && bx > cx) || (ax < bx && bx < cx))) {
      PASSMSG("minimum bracketed (func4)");
    } else {
      FAILMSG("minimum NOT bracketed (func4)");
    }
  }
  {
    std::cout << "\nExamine func3(x) = (x+1)*x*(x-1), for x in [0.0, 0.5]\n";
    ax = 0.0;
    bx = 0.5;
    rtt_min::mnbrak(ax, bx, cx, fa, fb, fc, func3);
    if ((fb < fa && fb < fc) && ((ax < bx && cx > bx) || (cx < bx && ax < bx))) {
      // minimum should be at x=0.577;
      PASSMSG("minimum bracketed (func3)");
    } else {
      FAILMSG("minimum NOT bracketed (func3)");
    }
  }
  { // Test the branch "RLimit parabolic u to maximum allowed value"
    std::cout << "\nExamine func5(x), for x in [-1.0, 0.5]\n";
    ax = -1.0;
    bx = 0.5;
    rtt_min::mnbrak(ax, bx, cx, fa, fb, fc, func5);
    if ((ax < bx && cx > bx) && (fa > fb && fc > fb)) {
      PASSMSG("minimum bracketed (func5)");
    } else {
      FAILMSG("minimum NOT bracketed (func5)");
    }
  }
  { // Test the branch "Reject parabolic u; use default magnification"
    std::cout << "\nExamine func6(x) = sin(x) +x/10, for x in [0.0, 2.0]\n";
    ax = 0.0;
    bx = 2.0;
    rtt_min::mnbrak(ax, bx, cx, fa, fb, fc, func6);
    if ((ax > cx && bx > cx) && rtt_dsxx::soft_equiv(fa, 0.0)) {
      PASSMSG("minimum bracketed (func6)");
    } else {
      FAILMSG("minimum NOT bracketed (func6)");
    }
  }
}

} // namespace

//------------------------------------------------------------------------------------------------//

int main(int argc, char *argv[]) {
  rtt_dsxx::ScalarUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    tstmnbrak(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstmnbrak.cc
//------------------------------------------------------------------------------------------------//
