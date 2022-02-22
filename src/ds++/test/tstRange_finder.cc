//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   ds++/test/tstRange_finder.cc
 * \author Mike Buksas
 * \date   Thu Feb  6 12:43:22 2003
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

#include "ds++/Range_Finder.hh"
#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"
#include <array>

using namespace std;
using namespace rtt_dsxx;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

void test_range_finder_left(UnitTest &ut) {

  array<double, 10> v = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0};

  vector<double> values(v.begin(), v.end());

  int index = Range_finder_left(v.begin(), v.end(), 1.5);
  FAIL_IF_NOT(index == 1);

  index = Range_finder_left(values.begin(), values.end(), 2.5);
  FAIL_IF_NOT(index == 2);

  // Check for equality at all values:
  for (int i = 0; i < 10; ++i) {
    index = Range_finder_left(v.begin(), v.end(), static_cast<double>(i));
    FAIL_IF_NOT(index == i);
  }

  // For equality with the last value, we should get n-1 with end catching:
  index = Range_finder_left_catch_end(v.begin(), v.end(), 9.0);
  FAIL_IF_NOT(index == 8);

  index = Range_finder_catch_end(v.begin(), v.end(), 9.0, RANGE_DIRECTION::LEFT);
  FAIL_IF_NOT(index == 8);

  //     index = Range_finder_left(v,v+10, 42.69);
  //     if (index != -1) ut.failure("test FAILS");

  //     index = Range_finder_left(v+5,v+10, 1.0);
  //     if (index != -1) ut.failure("test FAILS");

  array<double, 10> rv = {9.0, 8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0, 0.0};

  vector<double> rvalues(rv.begin(), rv.end());

  index = Range_finder(rvalues.rbegin(), rvalues.rend(), 5.5, RANGE_DIRECTION::LEFT);
  FAIL_IF_NOT(index == 5);

  index = Range_finder_left(rvalues.rbegin(), rvalues.rend(), 5.0);
  FAIL_IF_NOT(index == 5);

  //     index = Range_finder_left(rvalues.rbegin(), rvalues.rend(), 10.12);
  //     if (index != -1) ut.failure("test FAILS");
  using iterator = array<double, 10>::iterator;
  if (validate(pair<iterator, iterator>(rv.begin(), rv.begin()), rv.begin(), rv.end())) {
    ut.failure("validate FAILED to catch out of range result");
  } else {
    ut.passes("validate caught out of range result");
  }
  if (validate(pair<iterator, iterator>(rv.end(), rv.end()), rv.begin(), rv.end())) {
    ut.failure("validate FAILED to catch out of range result");
  } else {
    ut.passes("validate caught out of range result");
  }
  return;
}

//------------------------------------------------------------------------------------------------//
void test_range_finder_right(UnitTest &ut) {

  array<double, 10> v = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0};

  int index;

  // Check for equality at all values. Note that 0 comes back as interval -1
  // (e.g. out of range).
  for (int i = 1; i < 10; ++i) {
    index = Range_finder(v.begin(), v.end(), static_cast<double>(i), RANGE_DIRECTION::RIGHT);
    FAIL_IF_NOT(index == i - 1);
  }

  index = Range_finder_right_catch_end(v.begin(), v.end(), 0.0);
  FAIL_IF_NOT(index == 0);

  index = Range_finder_catch_end(v.begin(), v.end(), 0.0, RANGE_DIRECTION::RIGHT);
  FAIL_IF_NOT(index == 0);
  return;
}

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  ScalarUnitTest ut(argc, argv, release);
  try {
    test_range_finder_left(ut);
    test_range_finder_right(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstRange_finder.cc
//------------------------------------------------------------------------------------------------//
