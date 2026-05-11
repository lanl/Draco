//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   diagnostics/test/tstruntime_error.cc
 * \author Kent Budge
 * \brief  Test runtime_error function.
 * \note   Copyright (C) 2018-2026 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "dsxx/Release.hh"
#include "dsxx/UnitTest.hh"
#include "c4/C4_Functions.hh"
#include "c4/ParallelUnitTest.hh"
#include "diagnostics/runtime_check.hh"
#include <stdexcept>

namespace {

using namespace std;
using namespace rtt_dsxx;
using namespace rtt_c4;
using namespace rtt_diagnostics;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

void tstruntime_error(UnitTest &ut) {
  // Test passes
  try {
    rtt_diagnostics::runtime_check(true, "condition successful");
    PASSMSG("check true on all processors");
  } catch (...) {
    FAILMSG("check true on all processors");
  }

  try {
    runtime_check(false, "condition fails on all");
    FAILMSG("check false on all processors");
  } catch (std::runtime_error & /*error*/) {
    PASSMSG("check false on all processors throws on all as it should.");
  }

  try {
    runtime_check(rtt_c4::node() != 0, "condition fails on one processor");
    FAILMSG("check false on one processors");
  } catch (std::runtime_error & /*error*/) {
    PASSMSG("check false on one processor throws on all as it should.");
  }
}

} // namespace

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  rtt_c4::ParallelUnitTest ut(argc, argv, release);
  try {
    tstruntime_error(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstruntime_error.cc
//------------------------------------------------------------------------------------------------//
