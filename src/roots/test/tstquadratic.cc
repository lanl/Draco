//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   roots/test/tstquadratic.cc
 * \author Kent Budge
 * \date   Tue Jul  6 10:00:38 2004
 * \note   Copyright (C) 2010-2024 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "dsxx/Release.hh"
#include "dsxx/ScalarUnitTest.hh"
#include "dsxx/Soft_Equivalence.hh"
#include "roots/quadratic.hh"

using namespace std;
using namespace rtt_roots;
using namespace rtt_dsxx;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

int main(int argc, char *argv[]) {
  ScalarUnitTest ut(argc, argv, release);
  try {
    double r1, r2;
    quadratic(2.0, -2.0, -12.0, r1, r2);

    if (r1 > r2)
      swap(r1, r2);

    if (soft_equiv(r1, -2.0))
      PASSMSG("quadratic solve returned correct first root");
    else
      FAILMSG("quadratic solve returned INCORRECT first root");

    if (soft_equiv(r2, 3.0))
      PASSMSG("quadratic solve returned correct second root");
    else
      FAILMSG("quadratic solve returned INCORRECT second root");
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of testquadratic.cc
//------------------------------------------------------------------------------------------------//
