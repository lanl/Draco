//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   ds++/test/tstSoft_Equiv.cc
 * \author Thomas M. Evans
 * \date   Wed Nov  7 15:55:54 2001
 * \brief  Soft_Equiv header testing utilities.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"
#include "ds++/Soft_Equivalence.hh"
#include <array>
#include <deque>
#include <sstream>
#include <typeinfo>

using namespace std;
using rtt_dsxx::soft_equiv;
using rtt_dsxx::soft_equiv_deep;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

void test_soft_equiv_scalar(rtt_dsxx::ScalarUnitTest &ut) {

  double x = 0.9876543212345678;
  double y = 0.9876543212345678;

  FAIL_IF_NOT(soft_equiv(x, y, 1.e-16));
  FAIL_IF_NOT(soft_equiv(x, y));

  double z = 0.9876543212345679;
  FAIL_IF(soft_equiv(x, z, 1.e-16));

  double a = 0.987654321234;
  FAIL_IF_NOT(soft_equiv(x, a));

  a = 0.987654321233;
  FAIL_IF(soft_equiv(x, a));

  // checks for the new "reference=zero" coding 4aug00
  double zero = 0.0;
  FAIL_IF(soft_equiv(1.0e-10, zero));
  FAIL_IF(soft_equiv(-1.0e-10, zero));
  FAIL_IF_NOT(soft_equiv(-1.0e-35, zero));
  FAIL_IF_NOT(soft_equiv(1.0e-35, zero));

  return;
}

//------------------------------------------------------------------------------------------------//
void test_soft_equiv_container(rtt_dsxx::ScalarUnitTest &ut) {
  vector<double> values = {0.3247333291470, 0.3224333221471, 0.3324333522912};
  vector<double> const reference(values);

  if (soft_equiv(values.begin(), values.end(), reference.begin(), reference.end()))
    PASSMSG("Passed vector equivalence test.");
  else
    ITFAILS;

  // modify one value (delta < tolerance )
  values[1] += 1.0e-13;
  if (!soft_equiv(values.begin(), values.end(), reference.begin(), reference.end(), 1.e-13))
    PASSMSG("Passed vector equivalence precision test.");
  else
    ITFAILS;

  // Tests that compare 1D vector data to 1D array data.
  array<double, 3> v = {reference[0], reference[1], reference[2]};

  if (soft_equiv(v.begin(), v.end(), reference.begin(), reference.end()))
    PASSMSG("Passed vector-pointer equivalence test.");
  else
    ITFAILS;

  FAIL_IF_NOT(soft_equiv(reference.begin(), reference.end(), v.begin(), v.end()));

  // Check incompatible size
  FAIL_IF(soft_equiv(reference.begin(), reference.end(), v.begin() + 1, v.end()));

  // modify one value (delta < tolerance )
  v[1] += 1.0e-13;
  if (!soft_equiv(v.begin(), v.end(), reference.begin(), reference.end(), 1.e-13))
    PASSMSG("Passed vector-pointer equivalence precision test.");
  else
    ITFAILS;

  // C++ std::array containers
  std::array<double, 3> cppa_vals{{0.3247333291470, 0.3224333221471, 0.3324333522912}};
  if (soft_equiv(cppa_vals.begin(), cppa_vals.end(), reference.begin(), reference.end()))
    PASSMSG("Passed std::array<int,3> equivalence test.");
  else
    ITFAILS;

  // Try with a std::deque
  deque<double> d{};
  d.push_back(reference[0]);
  d.push_back(reference[1]);
  d.push_back(reference[2]);
  if (soft_equiv(d.begin(), d.end(), reference.begin(), reference.end()))
    PASSMSG("Passed deque<T> equivalence test.");
  else
    ITFAILS;

  return;
}

//------------------------------------------------------------------------------------------------//

void test_soft_equiv_deep_container(rtt_dsxx::ScalarUnitTest &ut) {

  vector<vector<double>> values = {{0.3247333291470, 0.3224333221471, 0.3324333522912},
                                   {0.3247333292470, 0.3224333222471, 0.3324333523912},
                                   {0.3247333293470, 0.3224333223471, 0.3324333524912}};
  vector<vector<double>> const reference = values;
  FAIL_IF_NOT(
      soft_equiv_deep<2>().equiv(values.begin(), values.end(), reference.begin(), reference.end()));

  // Soft_Equiv should still pass
  values[0][1] += 1.0e-13;
  FAIL_IF(soft_equiv_deep<2>().equiv(values.begin(), values.end(), reference.begin(),
                                     reference.end(), 1.e-13));

  // Compare C++ array to vector<vector<double>> data.
  // This cannot work because the C++ array is fundamentally a 1-D container.

  // double v[3];
  // v[0] = 0.3247333291470;
  // v[1] = 0.3224333221471;
  // v[2] = 0.3324333522912;
  // if (soft_equiv(&v[0], &v[3],
  //             reference.begin(), reference.end()))
  //     PASSMSG("Passed vector-pointer equivalence test.");
  // else
  //     ITFAILS;

  // if (!soft_equiv(reference.begin(), reference.end(), &v[0], &v[3]))
  //     ITFAILS;

  // Test 3-D array
  vector<vector<vector<double>>> const ref = {
      {{0.1, 0.2}, {0.3, 0.4}}, {{1.1, 1.2}, {1.3, 1.4}}, {{2.1, 2.2}, {2.3, 2.4}}};
  vector<vector<vector<double>>> val = ref;

  if (soft_equiv_deep<3>().equiv(val.begin(), val.end(), ref.begin(), ref.end()))
    PASSMSG("Passed vector<vector<vector<double>>> equivalence test.");
  else
    ITFAILS;

  if (!soft_equiv_deep<3>().equiv(val.begin(), val.end(), ref.begin() + 1, ref.end()))
    PASSMSG("Passed vector<vector<vector<double>>> equivalence test.");
  else
    ITFAILS;

  return;
}

//------------------------------------------------------------------------------------------------//
void test_vector_specialization(rtt_dsxx::ScalarUnitTest &ut) {
  double const epsilon(1.0e-27);
  {
    // 1-d vector comparison.
    std::vector<double> v(27, epsilon);
    std::vector<double> r(27, epsilon);
    FAIL_IF_NOT(soft_equiv(v, r));
  }
  {
    // 2-d vector comparison.
    std::vector<std::vector<double>> v(5);
    std::vector<std::vector<double>> r(5);
    for (size_t i = 0; i < 5; ++i) {
      v[i] = std::vector<double>(i * 2, epsilon);
      r[i] = std::vector<double>(i * 2, epsilon);
    }
    FAIL_IF_NOT(soft_equiv(v, r));
  }
  {
    // 3-d vector comparison.
    std::vector<std::vector<std::vector<double>>> v(5);
    std::vector<std::vector<std::vector<double>>> r(5);
    for (size_t i = 0; i < 5; ++i) {
      v[i].resize(3);
      r[i].resize(3);
      for (size_t j = 0; j < 3; ++j) {
        v[i][j] = std::vector<double>(j * 2, epsilon);
        r[i][j] = std::vector<double>(j * 2, epsilon);
      }
    }
    FAIL_IF_NOT(soft_equiv(v, r));
  }
  {
    // expect a failure for mismatched data.
    std::vector<double> v(27, 42.42);
    std::vector<double> r(27, 42.42);
    r[5] = 42.44; // mismatch value
    FAIL_IF(soft_equiv(v, r));
  }
  {
    // expect a failure for mismatched size.
    std::vector<double> v(27, 42.42);
    std::vector<double> r(7, 42.42);
    if (soft_equiv(v, r))
      FAILMSG("Comparing different size vectors should not pass!");
    else
      PASSMSG("Different size vectors are never equivalent.");
  }
  return;
}

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  rtt_dsxx::ScalarUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    // >>> UNIT TESTS
    test_soft_equiv_scalar(ut);
    test_soft_equiv_container(ut);
    test_soft_equiv_deep_container(ut);
    test_vector_specialization(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstSoft_Equiv.cc
//------------------------------------------------------------------------------------------------//
