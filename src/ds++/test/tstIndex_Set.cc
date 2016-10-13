//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ds++/test/tstIndex_Set.cc
 * \author Mike Buksas
 * \date   Thu Feb  2 13:46:36 2006
 * \brief  
 * \note   Copyright (C) 2016 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "ds++/Index_Set.hh"
#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"

using namespace std;
using namespace rtt_dsxx;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//
void test_index_set(rtt_dsxx::UnitTest &ut) {
  unsigned dimensions[] = {3, 4, 5};
  Index_Set<3, 1> box(dimensions);

  // Check the sizes and ranges of each dimensions
  if (box.get_size() != 60)
    ITFAILS;
  if (box.min_of_index() != 1)
    ITFAILS;
  if (box.max_of_index() != 60)
    ITFAILS;

  if (box.get_size(0) != 3)
    ITFAILS;
  if (box.min_of_index(0) != 1)
    ITFAILS;
  if (box.max_of_index(0) != 3)
    ITFAILS;

  if (box.get_size(1) != 4)
    ITFAILS;
  if (box.min_of_index(1) != 1)
    ITFAILS;
  if (box.max_of_index(1) != 4)
    ITFAILS;

  if (box.get_size(2) != 5)
    ITFAILS;
  if (box.min_of_index(2) != 1)
    ITFAILS;
  if (box.max_of_index(2) != 5)
    ITFAILS;

  if (box.limit_of_index(2, true) != 5)
    ITFAILS;
  if (box.limit_of_index(2, false) != 1)
    ITFAILS;
  if (box.limit_of_index(false) != 1)
    ITFAILS;

  // Test for indices in the total range and the range of each
  // dimension.
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

  unsigned indices[] = {4, 5, 6};
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

  if (box.get_size() != 1000)
    ITFAILS;

  dimensions[0] = 10;
  dimensions[1] = 2;
  dimensions[2] = 5;

  box.set_size(dimensions);

  if (box.get_size() != 100)
    ITFAILS;

  // Copy and comparison tests:
  // -------------------------

  // Make a copy
  Index_Set<3, 1> box2(box);

  // Test for equality and not-inequality
  if (box2 != box)
    ITFAILS;
  if (!(box2 == box))
    ITFAILS;

  // Resize the copy:
  box2.set_size(3);

  // Test for inequality and not-equality
  if (box2 == box)
    ITFAILS;
  if (!(box2 != box))
    ITFAILS;

  Index_Set<2, 0> csquare(4);

  // Check the sizes and ranges of each dimensions
  if (csquare.get_size() != 16)
    ITFAILS;
  if (csquare.min_of_index() != 0)
    ITFAILS;
  if (csquare.max_of_index() != 15)
    ITFAILS;

  if (ut.numFails == 0)
    PASSMSG("done with test_index_set().");
  return;
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[]) {
  rtt_dsxx::ScalarUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    test_index_set(ut);
  }
  UT_EPILOG(ut);
}

//---------------------------------------------------------------------------//
// end of tstIndex_Set.cc
//---------------------------------------------------------------------------//
