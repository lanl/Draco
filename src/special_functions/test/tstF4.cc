//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   special_functions/test/tstF4.cc
 * \author Kent Budge
 * \date   Tue Sep 21 11:57:47 2004
 * \note   Copyright (C) 2010-2024 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "dsxx/Release.hh"
#include "dsxx/ScalarUnitTest.hh"
#include "dsxx/Soft_Equivalence.hh"
#include "special_functions/F4.hh"
#include "units/PhysicalConstants.hh"
#include <fstream>

using namespace std;
using namespace rtt_dsxx;
using namespace rtt_sf;
using rtt_units::PI;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

void tstF4(UnitTest &ut) {
  double f3 = F4(-10.0);
  if (soft_equiv(f3, 4 * 3 * 2 * exp(-10.0) * (1 - exp(-10.0) / 32.0), 2e-6)) {
    ut.passes("correct F4 for -20.0");
  } else {
    ut.failure("NOT correct F4 for -20.0");
  }
  f3 = F4(1000.0);
  if (soft_equiv(f3, pow(1000.0, 5.0) / 5.0 + PI * PI * 4 * pow(1000.0, 3.0) / 6.0, 1.0e-10)) {
    ut.passes("correct F4 for 1000.0");
  } else {
    ut.failure("NOT correct F4 for 1000.0");
  }

  ofstream out("debug.dat");
  for (size_t j = 0; j < 300; ++j) {
    double eta = -10.0 + 0.1 * static_cast<double>(j);
    out << eta << ' ' << F4(eta) << endl;
  }
}

//------------------------------------------------------------------------------------------------//

int main(int argc, char *argv[]) {
  ScalarUnitTest ut(argc, argv, release);
  try {
    tstF4(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstF4.cc
//------------------------------------------------------------------------------------------------//
