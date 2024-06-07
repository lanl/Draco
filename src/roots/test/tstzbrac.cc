//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   roots/test/tstzbrac.cc
 * \author Kent Budge
 * \date   Tue Aug 17 15:24:48 2004
 * \brief  Test the zbrac function template
 * \note   Copyright (C) 2010-2024 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "dsxx/Release.hh"
#include "dsxx/ScalarUnitTest.hh"
#include "dsxx/fpe_trap.hh"
#include "roots/zbrac.hh"

using namespace std;
using namespace rtt_roots;
using namespace rtt_dsxx;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

double foo(double x) { return log(x); }

//------------------------------------------------------------------------------------------------//
//! Functor class representation of log function for testing zbrent and zbrac.
static class Zbrac_Test_Function {
public:
  double operator()(double x) const { return log(x); }
} zbrac_test_function;

//! Quadratic test
double cubic(double const x) {
  return -x * (x - 1) * (x - 2) / (-1 * (-1 - 1) * (-1 - 2)) +
         (x + 1) * (x - 1) * (x - 2) / ((0 + 1) * (0 - 1) * (0 - 2)) -
         (x + 1) * x * (x - 2) / ((1 + 1) * 1 * (1 - 2)) +
         (x + 1) * x * (x - 1) / ((2 + 1) * 2 * (2 - 1));
}

//------------------------------------------------------------------------------------------------//

void tstzbrac(UnitTest &ut) {
  double x1 = 0.1, x2 = 10.0;

  // Help out the compiler by defining function pointer signature.
  using fpdd = double (*)(double);
  fpdd log_fpdd = &foo;

  zbrac<fpdd>(log_fpdd, x1, x2);
  if (log(x1) * log(x2) < 0.0) {
    ut.passes("zbrac bracketed the zero of the log function");
  } else {
    ut.failure("zbrack did NOT bracket the zero of the log function");
  }

  // Try a slightly different starting search.

  x1 = 0.05, x2 = 0.1;
  zbrac<fpdd>(log_fpdd, x1, x2);
  if (log(x1) * log(x2) < 0.0) {
    ut.passes("zbrac bracketed the zero of the log function");
  } else {
    ut.failure("zbrack did NOT bracket the zero of the log function");
  }

  x1 = 9.0;
  x2 = 10.0;
  zbrac(zbrac_test_function, x1, x2);

  if (zbrac_test_function(x1) * zbrac_test_function(x2) < 0.0) {
    ut.passes("zbrac bracketed the zero of zbrac_test_function");
  } else {
    ut.failure("zbrack did NOT bracket the zero of zbrac_test_function");
  }

  // now do a function with more than one root
  x1 = -1;
  x2 = 1;
  zbrac<fpdd>(cubic, x1, x2);
  if (cubic(x1) * cubic(x2) < 0.0) {
    ut.passes("zbrac bracketed the zero of the cubic function");
  } else {
    ut.failure("zbrack did NOT bracket the zero of the cubic function");
  }
}

//------------------------------------------------------------------------------------------------//

int main(int argc, char *argv[]) {

  ScalarUnitTest ut(argc, argv, release);
  // This test includes a check for -NaN by design.  To avoid failure, do not run this test with
  // fpe_trap enabled.
  rtt_dsxx::fpe_trap fpeTrap(true);
  fpeTrap.disable();

  try {
    tstzbrac(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstzbrac.cc
//------------------------------------------------------------------------------------------------//
