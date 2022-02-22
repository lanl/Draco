//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   ds++/test/tstHomogeneous_New.cc
 * \author Kent Budge
 * \date   Tue Nov 28 09:17:23 2006
 * \brief  test the Homogeneous_New allocator class.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

#include "ds++/Homogeneous_New.hh"
#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"

using namespace std;
using namespace rtt_dsxx;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

void tstHomogeneous_New(UnitTest &ut) {
  { Homogeneous_New allocator(34); }
  ut.passes("Construction/destruction did not throw an exception");

  {
    Homogeneous_New allocator(56);
    void *first_object = allocator.allocate();
    cout << "First address: " << first_object << endl;
    void *second_object = allocator.allocate();
    cout << "Second address: " << second_object << endl;
    allocator.deallocate(first_object);
    first_object = allocator.allocate();
    cout << "Reallocated first address: " << first_object << endl;
    void *third_object = allocator.allocate();
    cout << "Third address: " << third_object << endl;
  }

  {
    Homogeneous_New allocator(29);
    allocator.reserve(7);
  }
}

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
