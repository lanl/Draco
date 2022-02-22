//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   ds++/test/tstscalar_unit_test.cc
 * \author Kent Grimmett Budge
 * \date   Tue Nov  6 13:19:40 2018
 * \brief  Test the test function for scalar case.
 * \note   Copyright (C) 2018-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"

using namespace rtt_dsxx;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

void test1(UnitTest &ut) { ut.passes("function test"); }

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  return do_scalar_unit_test(argc, argv, release, test1,
                             [](UnitTest &ut) { ut.passes("lambda test"); });
}

//------------------------------------------------------------------------------------------------//
// end of tstscalar_unit_test.cc
//------------------------------------------------------------------------------------------------//
