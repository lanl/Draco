//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   roots/test/tstzbrent.cc
 * \author Kent Budge
 * \brief  Test the zbrent function template
 * \note   Copyright (C) 2010-2026 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "dsxx/Release.hh"
#include "dsxx/ScalarUnitTest.hh"
#include "dsxx/Soft_Equivalence.hh"
#include "dsxx/UnitTest.hh"
#include "roots/zbrac.hh"
#include "roots/zbrent.hh"
#include <cmath> // log, exp
#include <exception>
#include <iostream>
#include <limits>

// cSpell:ignore fpdd, tstzbrent, zbrac, zbrent

namespace {

using namespace std;
using namespace rtt_dsxx;
using namespace rtt_roots;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//
double foo(double x) { return log(x); }
double foo_expm1(double x) { return exp(x) - 1.0; }

namespace {

//------------------------------------------------------------------------------------------------//
//! Functor class representation of log function for testing zbrent and zbrac.
class Zbrac_Test_Function {
public:
  double operator()(double x) const { return log(x); }
} zbrac_test_function;

} // namespace

//------------------------------------------------------------------------------------------------//
void tstzbrent(UnitTest &ut) {
  double x1 = 0.1, x2 = 10.0;

  // Help out the compiler by defining function pointer signature.
  using fpdd = double (*)(double);
  fpdd log_fpdd = &foo;
  zbrac<fpdd>(log_fpdd, x1, x2);

  double x_tol = numeric_limits<double>::epsilon();
  double f_tol = numeric_limits<double>::epsilon();

  if (soft_equiv(zbrent<fpdd>(log_fpdd, x1, x2, 100, x_tol, f_tol), 1.0)) {
    ut.passes("zbrent found the zero of the log function");
  } else {
    ut.failure("zbrent did NOT find the zero of the log function");
  }

  x1 = 0.1;
  x2 = 10.0;
  zbrac(zbrac_test_function, x1, x2);

  x_tol = numeric_limits<double>::epsilon();
  f_tol = numeric_limits<double>::epsilon();
  if (soft_equiv(zbrent(zbrac_test_function, x1, x2, 100, x_tol, f_tol), 1.0)) {
    ut.passes("zbrent found the zero of zbrac_test_function");
  } else {
    ut.failure("zbrent did NOT find the zero of zbrac_test_function");
  }

  x1 = 9.;
  x2 = 10.;

  fpdd expm1_fpdd = &foo_expm1;
  zbrac<fpdd>(expm1_fpdd, x1, x2);

  x_tol = numeric_limits<double>::epsilon();
  f_tol = numeric_limits<double>::epsilon();
  if (soft_equiv(zbrent<fpdd>(expm1_fpdd, x1, x2, 100, x_tol, f_tol), 0.0)) {
    ut.passes("zbrent found the zero of expm1");
  } else {
    ut.failure("zbrent did NOT find the zero of expm1");
  }
}

} // namespace

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  try {
    ScalarUnitTest ut(argc, argv, release);
    tstzbrent(ut);
  } catch (std::exception &err) {
    std::cout << "ERROR: While testing tstzbrent, " << err.what() << '\n';
    return 1;
  } catch (...) {
    std::cout << "ERROR: While testing tstzbrent, An unknown exception was thrown.\n";
    return 1;
  }
  return 0;
}

//------------------------------------------------------------------------------------------------//
// end of tstzbrent.cc
//------------------------------------------------------------------------------------------------//
