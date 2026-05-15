//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   dsxx/test/tstHomogeneous_New.cc
 * \author Kent Budge
 * \brief  test the Homogeneous_New allocator class.
 * \note   Copyright (C) 2024-2026 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "dsxx/Homogeneous_New.hh"
#include "dsxx/Release.hh"
#include "dsxx/ScalarUnitTest.hh"
#include "dsxx/UnitTest.hh"
#include <iostream>

namespace {

using namespace std;
using namespace rtt_dsxx;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

void tstHomogeneous_New(UnitTest &ut) {
  { Homogeneous_New const allocator(34); }
  ut.passes("Construction/destruction did not throw an exception");

  {
    Homogeneous_New allocator(56);
    void *first_object = allocator.allocate();
    cout << "First address: " << first_object << '\n';
    void *second_object = allocator.allocate();
    cout << "Second address: " << second_object << '\n';
    allocator.deallocate(first_object);
    first_object = allocator.allocate();
    cout << "Reallocated first address: " << first_object << '\n';
    void *third_object = allocator.allocate();
    cout << "Third address: " << third_object << '\n';
  }

  {
    Homogeneous_New allocator(29);
    allocator.reserve(7);
  }
}

} // namespace

//------------------------------------------------------------------------------------------------//

int main(int argc, char *argv[]) {
  ScalarUnitTest ut(argc, argv, release);
  try {
    tstHomogeneous_New(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstHomogeneous_New.cc
//------------------------------------------------------------------------------------------------//
