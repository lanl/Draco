//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   special_functions/test/tstFM12.cc
 * \author Kent Budge
 * \date   Tue Sep 21 11:57:47 2004
 * \note   Copyright (C) 2016-2020 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"
#include "ds++/Soft_Equivalence.hh"
#include "special_functions/FM12.hh"
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

void tstFM12(UnitTest &ut) {
  double f1 = FM12(-10.0);
  if (soft_equiv(f1, exp(-10.0 + gsl_sf_lngamma(0.5)) * (1 - exp(-10.0) / sqrt(2.)), 2e-6)) {
    ut.passes("correct FM12 for -10.0");
  } else {
    ut.failure("NOT correct FM12 for -10.0");
  }
  f1 = FM12(1000.0);
  if (soft_equiv(f1, pow(1000.0, 0.5) / 0.5 - PI * PI * 0.5 * pow(1000.0, -1.5) / 6.0, 1.0e-10)) {
    ut.passes("correct FM12 for 1000.0");
  } else {
    ut.failure("NOT correct FM12 for 1000.0");
  }

  ofstream out("debug.dat");
  for (size_t j = 0; j < 300; ++j) {
    double eta = -10.0 + 0.1 * static_cast<double>(j);
    out << eta << ' ' << FM12(eta) << endl;
  }
}

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  ScalarUnitTest ut(argc, argv, release);
  try {
    tstFM12(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstFM12.cc
//------------------------------------------------------------------------------------------------//
