//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   ode/test/tstquad.cc
 * \author Kent Budge
 * \date   Tue Sep 21 11:45:44 2004
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"
#include "ds++/Soft_Equivalence.hh"
#include "ode/quad.hh"
#include "ode/rkqs.hh"

using namespace std;
using namespace rtt_dsxx;
using namespace rtt_ode;

//------------------------------------------------------------------------------------------------//
// Helper function
//------------------------------------------------------------------------------------------------//
double foo_exp(double x) { return exp(x); }

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

void tstquad(UnitTest &ut) {
  using Rule = void (*)(vector<double> & y, const vector<double> &dydx, double &x,
                        const double htry, const double eps, const vector<double> &yscal,
                        double &hdid, double &hnext, Quad_To_ODE<double (*)(double)>);
  using fpdd = double (*)(double);
  fpdd exp_fpdd = &foo_exp;
  double eps = 1.0e-12;
  double const integral = rtt_ode::quad<fpdd, Rule>(exp_fpdd, 0.0, 1.0, eps, &rkqs);
  UT_MSG(soft_equiv(integral, exp(1.0) - 1.0, eps), "quad accurate");
}

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  ScalarUnitTest ut(argc, argv, release);
  try {
    tstquad(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstquad.cc
//------------------------------------------------------------------------------------------------//
