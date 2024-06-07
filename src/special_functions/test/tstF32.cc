//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   special_functions/test/tstF32.cc
 * \author Kent Budge
 * \date   Tue Sep 21 11:57:47 2004
 * \note   Copyright (C) 2010-2024 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "dsxx/Release.hh"
#include "dsxx/ScalarUnitTest.hh"
#include "dsxx/Soft_Equivalence.hh"
#include "special_functions/F32.hh"
#include "units/PhysicalConstants.hh"
#include <fstream>
#include <gsl/gsl_sf_gamma.h>

using namespace std;
using namespace rtt_dsxx;
using namespace rtt_sf;
using rtt_units::PI;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

void tstF32(UnitTest &ut) {
  double f1 = F32(-10.0);
  if (soft_equiv(f1, exp(-10.0 + gsl_sf_lngamma(2.5)) * (1 - exp(-10.0) / (4 * sqrt(2.))), 2e-6)) {
    ut.passes("correct F32 for -10.0");
  } else {
    ut.failure("NOT correct F32 for -10.0");
  }
  f1 = F32(1000.0);
  if (soft_equiv(f1, pow(1000.0, 2.5) / 2.5 + PI * PI * 1.5 * pow(1000.0, 0.5) / 6.0, 1.0e-10)) {
    ut.passes("correct F32 for 1000.0");
  } else {
    ut.failure("NOT correct F32 for 1000.0");
  }

  ofstream out("debug.dat");
  for (size_t j = 0; j < 300; ++j) {
    double eta = -10.0 + 0.1 * static_cast<double>(j);
    out << eta << ' ' << F32(eta) << endl;
  }
}

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  ScalarUnitTest ut(argc, argv, release);
  try {
    tstF32(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstF32.cc
//------------------------------------------------------------------------------------------------//
