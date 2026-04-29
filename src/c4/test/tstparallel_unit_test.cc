//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   c4/test/tstparallel_unit_test.cc
 * \author Kent Grimmett Budge
 * \brief  Test the test function for parallel case.
 * \note   Copyright (C) 2018-2026 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "dsxx/Release.hh"
#include "dsxx/UnitTest.hh"
#include "c4/ParallelUnitTest.hh"

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  return rtt_c4::do_parallel_unit_test(argc, argv, rtt_dsxx::release,
                                       [](rtt_dsxx::UnitTest &ut) { ut.passes("basic run"); });
}

//------------------------------------------------------------------------------------------------//
// end of tstparallel_unit_test.cc
//------------------------------------------------------------------------------------------------//
