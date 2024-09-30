//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   special_functions/test/tstF.cc
 * \author Kent Budge
 * \date   Tue Sep 21 11:57:47 2004
 * \brief  Unit test for F1 function.
 * \note   Copyright (C) 2024 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "dsxx/Release.hh"
#include "dsxx/ScalarUnitTest.hh"
#include "dsxx/Soft_Equivalence.hh"
#include "dsxx/UnitTest.hh"
#include "special_functions/F1.hh"
#include "special_functions/F12.hh"
#include "special_functions/F12inv.hh"
#include "special_functions/F2.hh"
#include "special_functions/F2inv.hh"
#include "special_functions/F3.hh"
#include "special_functions/F32.hh"
#include "special_functions/F4.hh"
#include "special_functions/FM12.hh"
#include "units/MathConstants.hh"
#include <cmath>   // exp, pow
#include <cstddef> // size_t
#include <fstream>
#include <gsl/gsl_sf_gamma.h>
#include <limits>

using namespace std;
using namespace rtt_dsxx;
using namespace rtt_sf;
using rtt_units::PI;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

void tstF1(UnitTest &ut) {
  double f1 = F1(-10.0);
  if (soft_equiv(f1, 1 * exp(-10.0) * (1 - exp(-10.0) / 4.0), 2e-6)) {
    PASSMSG("correct F1 for -10.0");
  } else {
    FAILMSG("NOT correct F1 for -10.0");
  }
  f1 = F1(1000.0);
  if (soft_equiv(f1, pow(1000.0, 2.0) / 2.0 + PI * PI * 1 * pow(1000.0, 0.0) / 6.0, 1.0e-10)) {
    PASSMSG("correct F1 for 1000.0");
  } else {
    FAILMSG("NOT correct F1 for 1000.0");
  }

  ofstream out("debug-F1.dat");
  for (size_t j = 0; j < 300; ++j) {
    double const eta = -10.0 + 0.1 * static_cast<double>(j);
    out << eta << ' ' << F1(eta) << '\n';
  }
}

//------------------------------------------------------------------------------------------------//
void tstF12(UnitTest &ut) {
  double f1 = F12(-10.0);
  if (soft_equiv(f1, exp(-10.0 + gsl_sf_lngamma(1.5)) * (1 - exp(-10.0) / (2 * sqrt(2.))), 2e-6)) {
    ut.passes("correct F12 for -10.0");
  } else {
    ut.failure("NOT correct F12 for -10.0");
  }
  f1 = F12(1000.0);
  if (soft_equiv(f1, pow(1000.0, 1.5) / 1.5 + PI * PI * 0.5 * pow(1000.0, -0.5) / 6.0, 1.0e-10)) {
    ut.passes("correct F12 for 1000.0");
  } else {
    ut.failure("NOT correct F12 for 1000.0");
  }

  ofstream out("debug-F12.dat");
  for (size_t j = 0; j < 300; ++j) {
    double const eta = -10.0 + 0.1 * static_cast<double>(j);
    out << eta << ' ' << F12(eta) << '\n';
  }
}

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

  ofstream out("debug-F32.dat");
  for (size_t j = 0; j < 300; ++j) {
    double const eta = -10.0 + 0.1 * static_cast<double>(j);
    out << eta << ' ' << F32(eta) << '\n';
  }
}

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

  ofstream out("debug-F4.dat");
  for (size_t j = 0; j < 300; ++j) {
    double const eta = -10.0 + 0.1 * static_cast<double>(j);
    out << eta << ' ' << F4(eta) << '\n';
  }
}

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

  ofstream out("debug-FM12.dat");
  for (size_t j = 0; j < 300; ++j) {
    double const eta = -10.0 + 0.1 * static_cast<double>(j);
    out << eta << ' ' << FM12(eta) << '\n';
  }
}

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  ScalarUnitTest ut(argc, argv, release);
  try {
    tstF1(ut);
    tstF12(ut);
    tstF3(ut);
    tstF32(ut);
    tstF4(ut);

    tstF12inv(ut);
    tstF2inv(ut);

    tstFM12(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstF.cc
//------------------------------------------------------------------------------------------------//
