//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   special_functions/test/tstF2inv.cc
 * \author Kent Budge
 * \date   Tue Sep 21 11:57:47 2004
 * \note   Copyright (C) 2010-2024 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "dsxx/Release.hh"
#include "dsxx/ScalarUnitTest.hh"
#include "dsxx/Soft_Equivalence.hh"
#include "special_functions/F2.hh"
#include "special_functions/F2inv.hh"
#include <limits>

using namespace std;
using namespace rtt_dsxx;
using namespace rtt_sf;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

void tstF2inv(UnitTest &ut) {
  // Nondegenerate limit
  double f12inv = F2inv(0.41223072e-8);
  if (soft_equiv(f12inv, -20.0, 2.0e-6)) {
    ut.passes("correct F2inv in nondegenerate limit");
  } else {
    ut.failure("NOT correct F2inv in nondegenerate limit");
  }
  f12inv = F2inv(0.90799344e-4);
  if (soft_equiv(f12inv, -10.0, 2.0e-6)) {
    ut.passes("correct F2inv in nondegenerate limit");
  } else {
    ut.failure("NOT correct F2inv in nondegenerate limit");
  }
  f12inv = F2inv(0.036551);
  if (soft_equiv(f12inv, -4.0, 2.0e-6)) {
    ut.passes("correct F2inv in nondegenerate limit");
  } else {
    ut.failure("NOT correct F2inv in nondegenerate limit");
  }
  f12inv = F2inv(2.821225);
  if (soft_equiv(f12inv, 0.5, 3.0e-6)) {
    ut.passes("correct F2inv in nondegenerate limit");
  } else {
    ut.failure("NOT correct F2inv in nondegenerate limit");
  }
  f12inv = F2inv(4.328723);
  if (soft_equiv(f12inv, 1.0, 5.0e-5)) {
    ut.passes("correct F2inv in degenerate limit");
  } else {
    ut.failure("NOT correct F2inv in degenerate limit");
  }
  f12inv = F2inv(165.31509);
  if (soft_equiv(f12inv, 7.5, 5.0e-6)) {
    ut.passes("correct F2inv in degenerate limit");
  } else {
    ut.failure("NOT correct F2inv in degenerate limit");
  }
  f12inv = F2inv(960.81184);
  if (soft_equiv(f12inv, 14.0, 5.0e-6)) {
    ut.passes("correct F2inv in degenerate limit");
  } else {
    ut.failure("NOT correct F2inv in degenerate limit");
  }
  f12inv = F2inv(2732.71153);
  if (soft_equiv(f12inv, 20.0, 5.0e-6)) {
    ut.passes("correct F2inv in degenerate limit");
  } else {
    ut.failure("NOT correct F2inv in degenerate limit");
  }
  f12inv = F2inv(9098.696);
  if (soft_equiv(f12inv, 30.0, 5.0e-6)) {
    ut.passes("correct F2inv in degenerate limit");
  } else {
    ut.failure("NOT correct F2inv in degenerate limit");
  }

  // Nondegenerate limit
  double f2 = F2(-20.0);
  if (soft_equiv(f2, 0.41223072e-8, 2.0e-6)) {
    ut.passes("correct F2 in nondegenerate limit");
  } else {
    ut.failure("NOT correct F2 in nondegenerate limit");
  }
  f2 = F2(-10.0);
  if (soft_equiv(f2, 0.90799344e-4, 2.0e-6)) {
    ut.passes("correct F2 in nondegenerate limit");
  } else {
    ut.failure("NOT correct F2 in nondegenerate limit");
  }
  f2 = F2(-4.0);
  if (soft_equiv(f2, 0.0365479, 2.0e-6)) {
    ut.passes("correct F2 in nondegenerate limit");
  } else {
    ut.failure("NOT correct F2 in nondegenerate limit");
  }
  f2 = F2(0.5);
  if (soft_equiv(f2, 2.82097, 3.0e-6)) {
    ut.passes("correct F2inv in nondegenerate limit");
  } else {
    ut.failure("NOT correct F2inv in nondegenerate limit");
  }
  f12inv = F2inv(4.328723);
  if (soft_equiv(f12inv, 1.0, 5.0e-5)) {
    ut.passes("correct F2inv in degenerate limit");
  } else {
    ut.failure("NOT correct F2inv in degenerate limit");
  }
  f2 = F2(7.5);
  if (soft_equiv(f2, 165.30012, 5.0e-6)) {
    ut.passes("correct F2 in degenerate limit");
  } else {
    ut.failure("NOT correct F2 in degenerate limit");
  }
  f12inv = F2inv(960.81184);
  if (soft_equiv(f12inv, 14.0, 5.0e-6)) {
    ut.passes("correct F2inv in degenerate limit");
  } else {
    ut.failure("NOT correct F2inv in degenerate limit");
  }
  f12inv = F2inv(2732.71153);
  if (soft_equiv(f12inv, 20.0, 5.0e-6)) {
    ut.passes("correct F2inv in degenerate limit");
  } else {
    ut.failure("NOT correct F2inv in degenerate limit");
  }
  f12inv = F2inv(9098.696);
  if (soft_equiv(f12inv, 30.0, 5.0e-6)) {
    ut.passes("correct F2inv in degenerate limit");
  } else {
    ut.failure("NOT correct F2inv in degenerate limit");
  }

  double const h = sqrt(numeric_limits<double>::epsilon());
  double mu, dmudf;
  f2 = F2(0.7);
  F2inv(f2, mu, dmudf);
  double fm = F2inv(f2 - h * f2);
  double fp = F2inv(f2 + h * f2);
  if (soft_equiv(mu, 0.7, 2e-4)) {
    ut.passes("correct F2inv");
  } else {
    ut.failure("NOT correct F2inv");
  }
  if (soft_equiv(dmudf, (fp - fm) / (2 * h * f2), 1e-8)) {
    ut.passes("correct F2inv deriv");
  } else {
    ut.failure("NOT correct F2inv deriv");
  }

  f2 = F2(9.7);
  F2inv(f2, mu, dmudf);
  fm = F2inv(f2 - h * f2);
  fp = F2inv(f2 + h * f2);
  if (soft_equiv(mu, 9.7, 5.0e-5)) {
    ut.passes("correct F2inv");
  } else {
    ut.failure("NOT correct F2inv");
  }
  if (soft_equiv(dmudf, (fp - fm) / (2 * h * f2), 5e-8)) {
    ut.passes("correct F2inv deriv");
  } else {
    ut.failure("NOT correct F2inv deriv");
  }
}

//------------------------------------------------------------------------------------------------//

int main(int argc, char *argv[]) {
  ScalarUnitTest ut(argc, argv, release);
  try {
    tstF2inv(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstF2inv.cc
//------------------------------------------------------------------------------------------------//
