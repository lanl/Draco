//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ds++/test/tstdbc.cc
 * \author Kent G. Budge
 * \date   Feb 18 2003
 * \brief  Copyright (C) 2016-2017 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "ds++/DracoMath.hh"
#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"
#include "ds++/dbc.hh"

using namespace std;
using namespace rtt_dsxx;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

class sum_predicate_Test_Predicate {
public:
  typedef double Return_Type;

  double operator()(const std::pair<double, char *> &p) { return p.first; }
};

//---------------------------------------------------------------------------//

void dbc_test(UnitTest &ut) {
  using std::pair;
  using std::vector;
  using rtt_dsxx::dim;
  using namespace std;

  if (!(abs(5.4) == 5.4) || !(abs(-2.1) == 2.1) || !(abs(0.0) == 0.0))
    ut.failure("abs function template FAILED");
  else
    ut.passes("abs function template ok");

  if (!(dim(2, 7) == 0) || !(dim(5, -3) == 8) || !(dim(4, 4) == 0))
    ut.failure("dim function template FAILED");
  else
    ut.passes("dim function template ok");

  double sum_test_array[6] = {1., 4., 3., 2., 5., 6.};

  if (is_monotonic_increasing(sum_test_array, sum_test_array + 6) ||
      !is_monotonic_increasing(sum_test_array, sum_test_array + 2))
    ut.failure("is_monotonic_increasing function template FAILED");
  else
    ut.passes("is_monotonic_increasing function template ok");

  // Ensure that the is_monotonic_increasing() function will return true if
  // there is only one data point.

  if (!is_monotonic_increasing(sum_test_array, sum_test_array + 1))
    ut.failure(string("is_monotonic_increasing function template ") +
               string("incorrectly reported length=1 container ") +
               string("non-monotonic."));
  else
    ut.passes(string("is_monotonic_increasing function template worked for ") +
              string("length=1 test."));

  if (is_strict_monotonic_increasing(sum_test_array, sum_test_array + 6) ||
      !is_strict_monotonic_increasing(sum_test_array, sum_test_array + 2))
    ut.failure("is_strict_monotonic_increasing function template FAILED");
  else
    ut.passes("is_strict_monotonic_increasing function template ok");

  // Ensure that the is_strict_monotonic_increasing() function will return
  // true if there is only one data point.

  if (is_strict_monotonic_increasing(sum_test_array, sum_test_array + 1))
    ut.passes(
        string("is_strict_monotonic_increasing function template worked ") +
        string("for length=1 test."));
  else
    ut.failure(
        string("is_strict_monotonic_increasing function template ") +
        string("incorrectly reported length=1 container non-monotonic."));

  if (is_strict_monotonic_decreasing(sum_test_array + 1, sum_test_array + 3) &&
      !is_strict_monotonic_decreasing(sum_test_array, sum_test_array + 6))
    ut.passes("is_strict_monotonic_decreasing function template ok");
  else
    ut.failure("is_strict_monotonic_decreasing function template FAILED");

  // Ensure that the is_strict_monotonic_decreasing() function will return
  // true if there is only one data point.

  if (is_strict_monotonic_decreasing(sum_test_array, sum_test_array + 1))
    ut.passes(string("is_strict_monotonic_decreasing function template ") +
              string("worked for length=1 test."));
  else
    ut.failure(string("is_strict_monotonic_decreasing function template ") +
               string("incorrectly reported length=1 container monotonic."));

  if (std::find_if(sum_test_array, sum_test_array + 6,
                   bind2nd(greater<double>(), 2.)) != sum_test_array + 1)
    ut.failure("std::bind2nd or std::greater function templates FAILED");
  else
    ut.passes("std::bind2nd or std::greater function templates ok");

// The preceeding tests whether binders are present in the C++
// implementation. The binders in <functional> supercede the old
// rtt_utils::exceeds predicate.

// Test badly formed numbers.
// Note: VS2013 with Nov 2013 CTP will not compile an explicit division by zero.
#ifndef WIN32
  if (!ut.fpe_trap_active) {
    double const Zero = 0.0;
    double const Infinity = 1.0 / Zero;
    if (rtt_dsxx::isInf(Infinity))
      ut.passes("isInfinity works on this platform");
    else
      ut.failure("isInfinity is problematic on this platform.");
    // This is optimized to double Nan(0.0) by Intel.
    // double Nan = Zero*Infinity;
    double const Nan = std::sqrt(-1.0);
    if (rtt_dsxx::isNan(Nan))
      ut.passes("isNaN works on this platform");
    else
      ut.failure("isNaN is problematic on this platform.");
  }
#endif

  // Check on symmetricity of matrix.
  vector<double> A(2 * 2);
  A[0 + 2 * 0] = 2.5;
  A[0 + 2 * 1] = 3.8;
  A[1 + 2 * 0] = 4.5;
  A[1 + 2 * 1] = 3.3;
  if (!is_symmetric_matrix(A, 2))
    ut.passes("detected nonsymmetric matrix");
  else
    ut.failure("did NOT detect nonsymmetric matrix");
  A[1 + 2 * 0] = A[0 + 2 * 1];
  if (is_symmetric_matrix(A, 2))
    ut.passes("passed symmetric matrix");
  else
    ut.failure("did NOT pass symmetric matrix");
  return;
}

//---------------------------------------------------------------------------//

void isFinite_test(UnitTest &ut) {
  double x(15.0);
  if (rtt_dsxx::isFinite(x))
    ut.passes("Correctly found x to be finite.");
  else
    ut.failure("Failed to find x to be finite.");
  return;
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[]) {
  ScalarUnitTest ut(argc, argv, release);
  try {
    dbc_test(ut);
    isFinite_test(ut);
  }
  UT_EPILOG(ut);
}

//---------------------------------------------------------------------------//
// end of tstdbc.cc
//---------------------------------------------------------------------------//
