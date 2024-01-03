//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   special_functions/test/tstF3.cc
 * \author Kent Budge
 * \date   Tue Sep 21 11:57:47 2004
 * \note   Copyright (C) 2010-2024 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//
#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"
#include "ds++/Soft_Equivalence.hh"
#include "special_functions/F3.hh"
#include "units/PhysicalConstants.hh"

using namespace std;
using namespace rtt_dsxx;
using namespace rtt_sf;
using rtt_units::PI;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

void tstF3(UnitTest &ut) {
  double f3 = F3(-10.0);
  if (soft_equiv(f3, 3 * 2 * exp(-10.0) * (1 - exp(-10.0) / 16.0), 2e-10)) {
    ut.passes("correct F3 for -20.0");
  } else {
    ut.failure("NOT correct F3 for -20.0");
  }
  f3 = F3(1000.0);
  if (soft_equiv(f3, pow(1000.0, 4.0) / 4.0 + PI * PI * 3 * pow(1000.0, 2.0) / 6.0, 1.0e-10)) {
    ut.passes("correct F3 for 1000.0");
  } else {
    ut.failure("NOT correct F3 for 1000.0");
  }
}

//------------------------------------------------------------------------------------------------//

int main(int argc, char *argv[]) {
  ScalarUnitTest ut(argc, argv, release);
  try {
    tstF3(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstF3.cc
//------------------------------------------------------------------------------------------------//
