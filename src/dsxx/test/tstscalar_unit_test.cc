//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   dsxx/test/tstscalar_unit_test.cc
 * \author Kent Grimmett Budge
 * \brief  Test the test function for scalar case.
 * \note   Copyright (C) 2024-2026 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "dsxx/Release.hh"
#include "dsxx/ScalarUnitTest.hh"
#include "dsxx/UnitTest.hh"

using namespace rtt_dsxx;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

namespace {
void test1(UnitTest &ut) { ut.passes("function test"); }
} // namespace

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  return do_scalar_unit_test(argc, argv, release, test1,
                             [](UnitTest &ut) { ut.passes("lambda test"); });
}

//------------------------------------------------------------------------------------------------//
// end of tstscalar_unit_test.cc
//------------------------------------------------------------------------------------------------//
