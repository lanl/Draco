//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   c4/test/tstparallel_unit_test.cc
 * \author Kent Grimmett Budge
 * \date   Tue Nov  6 13:19:40 2018
 * \brief  Test the test function for parallel case.
 * \note   Copyright (C) 2018-2024 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "dsxx/Release.hh"
#include "c4/ParallelUnitTest.hh"

using namespace rtt_dsxx;
using namespace rtt_c4;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  return do_parallel_unit_test(argc, argv, release, [](UnitTest &ut) { ut.passes("basic run"); });
}

//------------------------------------------------------------------------------------------------//
// end of tstparallel_unit_test.cc
//------------------------------------------------------------------------------------------------//
