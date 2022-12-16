//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   ds++/test/tstIndex_Counter.cc
 * \author Mike Buksas
 * \date   Wed Feb  1 08:58:48 2006
 * \brief  Unit test for Index_Counter
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "ds++/Index_Converter.hh"
#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"

using namespace std;
using namespace rtt_dsxx;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

void test_index_counter(rtt_dsxx::UnitTest &ut) {
  array<unsigned, 3> dimensions = {3, 4, 5};
  Index_Converter<3, 1> box(dimensions.data());
  Index_Converter<3, 1>::Counter it = box.counter();

  FAIL_IF_NOT(it.get_index() == 1);
  FAIL_IF_NOT(it.get_index(0) == 1);
  FAIL_IF_NOT(it.get_index(1) == 1);
  FAIL_IF_NOT(it.get_index(2) == 1);
  FAIL_IF_NOT(it.is_in_range());
  FAIL_IF_NOT(it.get_indices()[0] == 1);
  FAIL_IF_NOT(it.get_indices()[1] == 1);
  FAIL_IF_NOT(it.get_indices()[2] == 1);
  vector<unsigned> it_copy(3);
  it.get_indices(it_copy.begin());
  FAIL_IF_NOT(it_copy[0] == 1);
  FAIL_IF_NOT(it_copy[1] == 1);
  FAIL_IF_NOT(it_copy[2] == 1);

  ++it;

  FAIL_IF_NOT(it.get_index() == 2);
  FAIL_IF_NOT(it.get_index(0) == 2);
  FAIL_IF_NOT(it.get_index(1) == 1);
  FAIL_IF_NOT(it.get_index(2) == 1);
  FAIL_IF_NOT(it.is_in_range());

  --it;
  --it;

  FAIL_IF(it.is_in_range());

  if (ut.numFails == 0)
    PASSMSG("done with test_index_counter().");
  return;
}

//------------------------------------------------------------------------------------------------//
void test_looping(rtt_dsxx::UnitTest &ut) {
  array<unsigned, 3> dimensions = {3, 4, 5};
  Index_Converter<3, 1> box(dimensions.data());

  int index = 1;
  for (Index_Counter<3, 1> it(box); it.is_in_range(); ++it) {
    const int it_index = it.get_index();

    // Check the returned index against a manual count.
    FAIL_IF_NOT(it_index == index++);

    // Check the first and last index directly.
    FAIL_IF_NOT((it_index - 1) % 3 + 1 == it.get_index(0));
    FAIL_IF_NOT((it_index - 1) / 12 + 1 == it.get_index(2));
  }
  return;
}

//------------------------------------------------------------------------------------------------//
void test_next_index(rtt_dsxx::UnitTest &ut) {
  array<unsigned, 3> dimensions = {3, 4, 5};
  Index_Converter<3, 1> box(dimensions.data());
  Index_Counter<3, 1> it = box.counter();

  FAIL_IF_NOT(it.get_index() == 1);
  FAIL_IF_NOT(box.get_next_index(it, 1) == -1);
  FAIL_IF_NOT(box.get_next_index(it, 2) == 2);
  FAIL_IF_NOT(box.get_next_index(it, 3) == -1);
  FAIL_IF_NOT(box.get_next_index(it, 4) == 4);
  FAIL_IF_NOT(box.get_next_index(it, 5) == -1);
  FAIL_IF_NOT(box.get_next_index(it, 6) == 13);
  return;
}

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  rtt_dsxx::ScalarUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    test_index_counter(ut);
    test_looping(ut);
    test_next_index(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstIndex_Counter.cc
//------------------------------------------------------------------------------------------------//
