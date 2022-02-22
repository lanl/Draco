//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   ds++/test/tstIndex_Set.cc
 * \author Mike Buksas
 * \date   Thu Feb  2 13:46:36 2006
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

#include "ds++/Index_Set.hh"
#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"

using namespace std;
using namespace rtt_dsxx;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//
void test_index_set(rtt_dsxx::UnitTest &ut) {
  array<unsigned, 3> dimensions = {3, 4, 5};
  Index_Set<3, 1> box(dimensions.data());

  // Check the sizes and ranges of each dimensions
  FAIL_IF_NOT(box.get_size() == 60);
  FAIL_IF_NOT(box.min_of_index() == 1);
  FAIL_IF_NOT(box.max_of_index() == 60);

  FAIL_IF_NOT(box.get_size(0) == 3);
  FAIL_IF_NOT(box.min_of_index(0) == 1);
  FAIL_IF_NOT(box.max_of_index(0) == 3);

  FAIL_IF_NOT(box.get_size(1) == 4);
  FAIL_IF_NOT(box.min_of_index(1) == 1);
  FAIL_IF_NOT(box.max_of_index(1) == 4);

  FAIL_IF_NOT(box.get_size(2) == 5);
  FAIL_IF_NOT(box.min_of_index(2) == 1);
  FAIL_IF_NOT(box.max_of_index(2) == 5);

  FAIL_IF_NOT(box.limit_of_index(2, true) == 5);
  FAIL_IF_NOT(box.limit_of_index(2, false) == 1);
  FAIL_IF_NOT(box.limit_of_index(false) == 1);

  // Test for indices in the total range and the range of each dimension.
  if (box.index_in_range(0))
    ITFAILS;
  if (!box.index_in_range(1))
    ITFAILS;
  if (!box.index_in_range(60))
    ITFAILS;
  if (box.index_in_range(61))
    ITFAILS;

  if (box.index_in_range(0, 0))
    ITFAILS;
  if (!box.index_in_range(1, 0))
    ITFAILS;
  if (!box.index_in_range(3, 0))
    ITFAILS;
  if (box.index_in_range(4, 0))
    ITFAILS;

  if (box.index_in_range(0, 0))
    ITFAILS;
  if (!box.index_in_range(1, 0))
    ITFAILS;
  if (!box.index_in_range(3, 0))
    ITFAILS;
  if (box.index_in_range(4, 0))
    ITFAILS;

  array<unsigned, 3> indices = {4, 5, 6};
  if (box.indices_in_range(indices))
    ITFAILS;

  // Test the functions for vetting direction and dimension arguments.
  if (box.direction_okay(0))
    ITFAILS;
  if (!box.direction_okay(1))
    ITFAILS;
  if (!box.direction_okay(6))
    ITFAILS;
  if (box.direction_okay(7))
    ITFAILS;

  // if (box.dimension_okay(-1)) ITFAILS;
  if (!box.dimension_okay(0))
    ITFAILS;
  if (!box.dimension_okay(2))
    ITFAILS;
  if (box.dimension_okay(3))
    ITFAILS;

  // Resize the object and repeat some of the tests:

  // Make a uniform array 10x10x10
  box.set_size(10);

  FAIL_IF_NOT(box.get_size() == 1000);

  dimensions[0] = 10;
  dimensions[1] = 2;
  dimensions[2] = 5;

  box.set_size(dimensions.data());

  FAIL_IF_NOT(box.get_size() == 100);

  // Copy and comparison tests:
  // -------------------------

  // Make a copy
  Index_Set<3, 1> box2(box);

  // Test for equality and not-inequality
  FAIL_IF_NOT(box2 == box);
  if (!(box2 == box))
    ITFAILS;

  // Resize the copy:
  box2.set_size(3);

  // Test for inequality and not-equality
  FAIL_IF(box2 == box);
  FAIL_IF_NOT(box2 != box);

  Index_Set<2, 0> csquare(4);

  // Check the sizes and ranges of each dimensions
  FAIL_IF_NOT(csquare.get_size() == 16);
  FAIL_IF_NOT(csquare.min_of_index() == 0);
  FAIL_IF_NOT(csquare.max_of_index() == 15);

  if (ut.numFails == 0)
    PASSMSG("done with test_index_set().");
  return;
}

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  rtt_dsxx::ScalarUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    test_index_set(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstIndex_Set.cc
//------------------------------------------------------------------------------------------------//
