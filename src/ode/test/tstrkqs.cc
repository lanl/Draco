//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   ode/test/tstrkqs.cc
 * \author Kent Budge
 * \brief  Test the rkqs integrator function template.
 * \note   Copyright (C) 2010-2026 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "dsxx/Release.hh"
#include "dsxx/ScalarUnitTest.hh"
#include "dsxx/UnitTest.hh"
#include "ode/rkqs.hh"
#include <cmath>
#include <limits>
#include <numbers>
#include <vector>

namespace {

using namespace std;
using namespace rtt_ode;
using namespace rtt_dsxx;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

// Instead of declaring ederivs 'static', use an unnamed namespace to limit the scope to this file.
void ederivs(double /*x*/, const vector<double> &y, vector<double> &dydx) { dydx[0] = y[0]; }

//------------------------------------------------------------------------------------------------//
void tstrkqs(UnitTest &ut) {
  vector<double> y(1, 1.0);
  vector<double> dydx(1, 1.0);
  double x = 0.0;
  vector<double> yout, yerr;
  rkck(y, dydx, x, 1.0, yout, yerr, ederivs);

  if (fabs(yout[0] - std::numbers::e) > 1.0e-2) {
    FAILMSG("rkck NOT accurate");
  } else {
    PASSMSG("rkck accurate");
  }

  y[0] = 1.0;
  dydx[0] = 1.0;
  x = 0.0;
  vector<double> const yscal(1, 1.0);
  double hdid, hnext;
  rkqs(y, dydx, x, 1.0, numeric_limits<double>::epsilon(), yscal, hdid, hnext, ederivs);

  if (fabs(y[0] - exp(x)) > 1.0e-9) {
    FAILMSG("rkqs NOT accurate");
  } else {
    PASSMSG("rkqs accurate");
  }
}

} // namespace

//------------------------------------------------------------------------------------------------//

int main(int argc, char *argv[]) {
  ScalarUnitTest ut(argc, argv, release);
  try {
    tstrkqs(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstrkqs.cc
//------------------------------------------------------------------------------------------------//
