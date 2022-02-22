//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   ode/test/tstrkqs.cc
 * \author Kent Budge
 * \date   Tue Sep 21 09:28:30 2004
 * \brief  Test the rkqs integrator function template.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

//------------------------------------------------------------------------------------------------//

#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"
#include "ode/rkqs.hh"

#include <cmath>
#include <limits>

using namespace std;
using namespace rtt_ode;
using namespace rtt_dsxx;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

static void ederivs(double /*x*/, const vector<double> &y, vector<double> &dydx) { dydx[0] = y[0]; }

void tstrkqs(UnitTest &ut) {
  vector<double> y(1, 1.0);
  vector<double> dydx(1, 1.0);
  double x = 0.0;
  vector<double> yout, yerr;
  rkck(y, dydx, x, 1.0, yout, yerr, ederivs);

  if (fabs(yout[0] - exp(1.0)) > 1.0e-2) {
    ut.failure("rkck NOT accurate");
  } else {
    ut.passes("rkck accurate");
  }

  y[0] = 1.0;
  dydx[0] = 1.0;
  x = 0.0;
  vector<double> yscal(1, 1.0);
  double hdid, hnext;
  rkqs(y, dydx, x, 1.0, numeric_limits<double>::epsilon(), yscal, hdid, hnext, ederivs);

  if (fabs(y[0] - exp(x)) > 1.0e-9) {
    ut.failure("rkqs NOT accurate");
  } else {
    ut.passes("rkqs accurate");
  }
}

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
