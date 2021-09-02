//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   ds++/test/tstdbc.cc
 * \author Kent G. Budge
 * \date   Feb 18 2003
 * \brief  Copyright (C) 2010-2021 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "ds++/DracoMath.hh"
#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"
#include "ds++/dbc.hh"
#include <array>

using namespace std;
using namespace rtt_dsxx;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

class sum_predicate_Test_Predicate {
public:
  using Return_Type = double;

  double operator()(const std::pair<double, char *> &p) { return p.first; }
};

//------------------------------------------------------------------------------------------------//
void dbc_test(UnitTest &ut) {
  using rtt_dsxx::dim;
  using std::pair;
  using std::vector;
  using namespace std;

  if (abs(abs(5.4) - 5.4) > std::numeric_limits<double>::epsilon() ||
      abs(abs(-2.1) - 2.1) > std::numeric_limits<double>::epsilon() ||
      abs(abs(0.0) - 0.0) > std::numeric_limits<double>::epsilon())
    FAILMSG("abs function template FAILED");
  else
    PASSMSG("abs function template ok");

  if (abs(dim(2, 7) - 0) > std::numeric_limits<int>::epsilon() ||
      abs(dim(5, -3) - 8) > std::numeric_limits<int>::epsilon() ||
      abs(dim(4, 4) - 0) > std::numeric_limits<int>::epsilon())
    FAILMSG("dim function template FAILED");
  else
    PASSMSG("dim function template ok");

  array<double, 6> sum_test_array = {1., 4., 3., 2., 5., 6.};

  if (is_monotonic_increasing(sum_test_array.begin(), sum_test_array.end()) ||
      !is_monotonic_increasing(sum_test_array.begin(), sum_test_array.begin() + 2))
    FAILMSG("is_monotonic_increasing function template FAILED");
  else
    PASSMSG("is_monotonic_increasing function template ok");

  // Ensure that the is_monotonic_increasing() function will return true if there is only one data
  // point.

  if (!is_monotonic_increasing(sum_test_array.begin(), sum_test_array.begin() + 1))
    FAILMSG(string("is_monotonic_increasing function template ") +
            "incorrectly reported length=1 container non-monotonic.");
  else
    PASSMSG(string("is_monotonic_increasing function template worked for ") +
            string("length=1 test."));

  if (is_strict_monotonic_increasing(sum_test_array.begin(), sum_test_array.end()) ||
      !is_strict_monotonic_increasing(sum_test_array.begin(), sum_test_array.begin() + 2))
    FAILMSG("is_strict_monotonic_increasing function template FAILED");
  else
    PASSMSG("is_strict_monotonic_increasing function template ok");

  // Ensure that the is_strict_monotonic_increasing() function will return true if there is only one
  // data point.

  if (is_strict_monotonic_increasing(sum_test_array.begin(), sum_test_array.begin() + 1))
    PASSMSG(string("is_strict_monotonic_increasing function template worked ") +
            string("for length=1 test."));
  else
    FAILMSG(string("is_strict_monotonic_increasing function template ") +
            string("incorrectly reported length=1 container non-monotonic."));

  if (is_strict_monotonic_decreasing(sum_test_array.begin() + 1, sum_test_array.begin() + 3) &&
      !is_strict_monotonic_decreasing(sum_test_array.begin(), sum_test_array.end()))
    PASSMSG("is_strict_monotonic_decreasing function template ok");
  else
    FAILMSG("is_strict_monotonic_decreasing function template FAILED");

  // Ensure that the is_strict_monotonic_decreasing() function will return true if there is only one
  // data point.

  if (is_strict_monotonic_decreasing(sum_test_array.begin(), sum_test_array.begin() + 1))
    PASSMSG(string("is_strict_monotonic_decreasing function template ") +
            string("worked for length=1 test."));
  else
    FAILMSG(string("is_strict_monotonic_decreasing function template ") +
            string("incorrectly reported length=1 container monotonic."));

  if (std::find_if(sum_test_array.begin(), sum_test_array.end(),
                   [](double x) { return x > 2.0; }) != sum_test_array.begin() + 1)
    FAILMSG("std::bind2nd or std::greater function templates FAILED");
  else
    PASSMSG("std::bind2nd or std::greater function templates ok");

// The preceeding tests whether binders are present in the C++ implementation. The binders in
// <functional> supercede the old rtt_utils::exceeds predicate.

// Test badly formed numbers.
// Note: VS2013 with Nov 2013 CTP will not compile an explicit division by zero.
#ifndef _MSC_VER
  if (!ut.fpe_trap_active) {
    double const Zero = 0.0;
    double const Infinity = 1.0 / Zero;
    if (rtt_dsxx::isInf(Infinity))
      PASSMSG("isInfinity works on this platform");
    else
      FAILMSG("isInfinity is problematic on this platform.");
    // This is optimized to double Nan(0.0) by Intel.
    // double Nan = Zero*Infinity;
    double const Nan = std::sqrt(-1.0);
    if (rtt_dsxx::isNan(Nan))
      PASSMSG("isNaN works on this platform");
    else
      FAILMSG("isNaN is problematic on this platform.");
  }
#endif

  // Check on symmetricity of matrix.
  vector<double> A(2 * 2);
  A[0 + 2 * 0] = 2.5;
  A[0 + 2 * 1] = 3.8;
  A[1 + 2 * 0] = 4.5;
  A[1 + 2 * 1] = 3.3;
  UT_MSG(!is_symmetric_matrix(A, 2), "detected nonsymmetric matrix");

  A[1 + 2 * 0] = A[0 + 2 * 1];
  UT_MSG(is_symmetric_matrix(A, 2), "passed symmetric matrix");

  return;
}

//------------------------------------------------------------------------------------------------//
void isFinite_test(UnitTest &ut) {
  UT_MSG(rtt_dsxx::isFinite(15.0), "Correctly found 15.0 to be finite.");
  return;
}

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  ScalarUnitTest ut(argc, argv, release);
  try {
    dbc_test(ut);
    isFinite_test(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstdbc.cc
//------------------------------------------------------------------------------------------------//
