//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   c4/test/tstGlobalTimer.cc
 * \author Kelly G. Thompson <kgt@lanl.gov>
 * \date   Thursday, Sep 13, 2018, 10:34 am
 * \brief  Test global timing functions in C4.
 * \note   Copyright (C) 2018-2022 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

#include "c4/Global_Timer.hh"
#include "c4/ParallelUnitTest.hh"
#include "ds++/Release.hh"
#include <chrono>
#include <sstream>
#include <thread>

//------------------------------------------------------------------------------------------------//
// TESTS
//
// Note that rtt_c4::Global_Timer::set_global_activity(bool) seems to conflict
// with rtt_c4::Timer, so this check was moved to an independent test.
//------------------------------------------------------------------------------------------------//

// This test fails to run.  Probably incorrect Use of Global_Timer.
void test_Global_Timer(rtt_dsxx::UnitTest &ut) {
  using namespace std::chrono_literals;

  std::cout << "Starting tstTime::test_Global_Timer tests..." << std::endl;

  rtt_c4::Global_Timer::set_global_activity(true);

  // sleep so the timer has something to report.
  std::this_thread::sleep_for(500ms);

  // This seems to report nothing, but it runs w/o error.
  rtt_c4::Global_Timer::report_all(std::cout);

  rtt_c4::Global_Timer::set_global_activity(false);

  {
    // Test dynamic object creation and destruction
    auto *ptr_gt = new rtt_c4::Global_Timer("foo");
    delete ptr_gt;
  }

  PASSMSG("done");
  return;
}

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  rtt_c4::ParallelUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    test_Global_Timer(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstTime.cc
//------------------------------------------------------------------------------------------------//
