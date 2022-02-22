//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   c4/test/tstCompare.cc
 * \author Mike Buksas
 * \date   Thu May  1 14:47:00 2008
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

#include "c4/Compare.hh"
#include "c4/ParallelUnitTest.hh"
#include "ds++/Release.hh"
#include <sstream>

using namespace std;
using namespace rtt_c4;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

template <typename T>
void test_equivalence(rtt_dsxx::UnitTest &ut, const T value, const T alt_value) {
  T local_value = value;

  // Test requires more than one node:
  if (rtt_c4::nodes() > 1) {
    // at this point all processors should have the same value
    FAIL_IF_NOT(check_global_equiv(local_value));

    // now change the first processor's value
    if (rtt_c4::node() == 0)
      local_value = alt_value;

    if (rtt_c4::node() > 0) {
      FAIL_IF_NOT(check_global_equiv(local_value));
    } else {
      FAIL_IF(check_global_equiv(local_value));
    }

    // Reset to given value
    local_value = value;
    FAIL_IF_NOT(check_global_equiv(local_value));

    // Change the last processor:
    if (rtt_c4::node() == rtt_c4::nodes() - 1)
      local_value = alt_value;

    if (rtt_c4::node() == rtt_c4::nodes() - 2) {
      FAIL_IF(check_global_equiv(local_value));
    } else {
      FAIL_IF_NOT(check_global_equiv(local_value));
    }
  }

  // Reset to given value
  local_value = value;
  FAIL_IF_NOT(check_global_equiv(local_value));

  // Test valid on two nodes or more:
  if (rtt_c4::nodes() > 2) {
    // Change a middle value
    if (rtt_c4::node() == rtt_c4::nodes() / 2)
      local_value = alt_value;

    if (rtt_c4::node() == rtt_c4::nodes() / 2 - 1) {
      FAIL_IF(check_global_equiv(local_value));
    } else if (rtt_c4::node() == rtt_c4::nodes() / 2) {
      FAIL_IF(check_global_equiv(local_value));
    } else {
      FAIL_IF_NOT(check_global_equiv(local_value));
    }
  }

  // Reset
  local_value = value;
  FAIL_IF_NOT(check_global_equiv(local_value));

  // Check 1 node. trivial, but check anyway.
  if (rtt_c4::nodes() == 1) {
    local_value = alt_value;
    FAIL_IF_NOT(check_global_equiv(local_value));
  }

  if (ut.numFails == 0) {
    std::ostringstream msg;
    msg << "No failures detected for test_equivalence(ut," << value << "," << alt_value << ").";
    PASSMSG(msg.str());
  }
  return;
}

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  rtt_c4::ParallelUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    // test global equivalences
    test_equivalence(ut, 10, 11);           // int
    test_equivalence(ut, 10.0001, 11.0001); // double
    test_equivalence(ut, 10.0001, 10.0002); // double
    test_equivalence(ut, static_cast<unsigned long long>(10000000000),
                     static_cast<unsigned long long>(200000000000));
    test_equivalence(ut, static_cast<long long>(10000000000), static_cast<long long>(200000000000));
    test_equivalence(ut, static_cast<long>(1000000), static_cast<long>(2000000));
    test_equivalence(ut, static_cast<unsigned long>(1000000), static_cast<unsigned long>(2000000));
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstCompare.cc
//------------------------------------------------------------------------------------------------//
