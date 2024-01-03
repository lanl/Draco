//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   special_functions/test/tstF12inv.cc
 * \author Kent Budge
 * \date   Tue Sep 21 11:57:47 2004
 * \note   Copyright (C) 2010-2024 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"
#include "ds++/Soft_Equivalence.hh"
#include "special_functions/F12.hh"
#include "special_functions/F12inv.hh"
#include <cmath>
#include <limits>

using namespace std;
using namespace rtt_dsxx;
using namespace rtt_sf;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

void tstF12inv(UnitTest &ut) {
  double f12 = F12(1.3);
  double f12inv = F12inv(f12);
  if (soft_equiv(f12inv, 1.3, 1.0e-9)) {
    ut.passes("correct F12inv");
  } else {
    ut.failure("NOT correct F12inv");
  }
  f12 = F12(9.7);
  f12inv = F12inv(f12);
  if (soft_equiv(f12inv, 9.7, 5.0e-9)) {
    ut.passes("correct F12inv");
  } else {
    ut.failure("NOT correct F12inv");
  }

  double const h = sqrt(numeric_limits<double>::epsilon());
  double mu, dmudf;
  f12 = F12(1.3);
  F12inv(f12, mu, dmudf);
  double fm = F12inv(f12 - h * f12);
  double fp = F12inv(f12 + h * f12);
  if (soft_equiv(mu, 1.3, 1.0e-9)) {
    ut.passes("correct F12inv");
  } else {
    ut.failure("NOT correct F12inv");
  }
  if (soft_equiv(dmudf, (fp - fm) / (2. * h * f12), 1.e-8)) {
    ut.passes("correct F12inv deriv");
  } else {
    ut.failure("NOT correct F12inv deriv");
  }

  f12 = F12(9.7);
  F12inv(f12, mu, dmudf);
  fm = F12inv(f12 - h * f12);
  fp = F12inv(f12 + h * f12);
  if (soft_equiv(mu, 9.7, 5.0e-9)) {
    ut.passes("correct F12inv");
  } else {
    ut.failure("NOT correct F12inv");
  }
  if (soft_equiv(dmudf, (fp - fm) / (2 * h * f12), 1e-8)) {
    ut.passes("correct F12inv deriv");
  } else {
    ut.failure("NOT correct F12inv deriv");
  }
}

//------------------------------------------------------------------------------------------------//

int main(int argc, char *argv[]) {
  ScalarUnitTest ut(argc, argv, release);
  try {
    tstF12inv(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstF12inv.cc
//------------------------------------------------------------------------------------------------//
