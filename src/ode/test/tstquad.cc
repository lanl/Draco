//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ode/test/tstquad.cc
 * \author Kent Budge
 * \date   Tue Sep 21 11:45:44 2004
 * \note   Copyright 2016-2017 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"
#include "ds++/Soft_Equivalence.hh"
#include "ode/quad.hh"
#include "ode/rkqs.hh"

using namespace std;
using namespace rtt_dsxx;
using namespace rtt_ode;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//
double foo_exp(double x) { return exp(x); }

void tstquad(UnitTest &ut) {
  typedef void (*Rule)(vector<double> & y, const vector<double> &dydx,
                       double &x, const double htry, const double eps,
                       const vector<double> &yscal, double &hdid, double &hnext,
                       Quad_To_ODE<double (*)(double)>);

  typedef double (*fpdd)(double);
  fpdd exp_fpdd = &foo_exp;

  double eps = 1.0e-12;
  double integral = rtt_ode::quad<fpdd, Rule>(exp_fpdd, 0.0, 1.0, eps, &rkqs);

  if (!soft_equiv(integral, exp(1.0) - 1.0, 1.0e-12)) {
    ut.failure("quad NOT accurate");
  } else {
    ut.passes("quad accurate");
  }
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[]) {
  ScalarUnitTest ut(argc, argv, release);
  try {
    tstquad(ut);
  }
  UT_EPILOG(ut);
}

//---------------------------------------------------------------------------//
// end of tstquad.cc
//---------------------------------------------------------------------------//
