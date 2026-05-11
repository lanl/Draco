//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   diagnostics/test/tstProcmon.cc
 * \author Kelly Thompson
 * \brief  Procmon test.
 * \note   Copyright (C) 2013-2026 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "dsxx/DracoStrings.hh"
#include "dsxx/Release.hh"
#include "dsxx/ScalarUnitTest.hh"
#include "dsxx/UnitTest.hh"
#include "diagnostics/Procmon.hh"
#include <cstdlib>
#include <iostream>
#include <numbers>
#include <sstream>
#include <string>
#include <vector>

// cSpell:ignore vmrss, vdbl

namespace {

using namespace std;
using namespace rtt_diagnostics;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

void tst_procmon_basic(rtt_dsxx::UnitTest &ut) {
  std::cout << "Running tst_procmon_basic(ut)...\n";

  int const mynode(0);
  double vmrss_orig(-1.0);
  double vmrss_bigger(-1.0);
  double vmrss_end(-1.0);

  // Run the function and ensure the output has some required fields.
  {
    std::ostringstream msg;
    procmon_resource_print("tst_procmon_use_01", mynode, msg);
    std::string const line = msg.str();
    std::vector<std::string> tokens = rtt_dsxx::tokenize(line, " \t");

    std::cout << msg.str() << '\n';
    for (size_t i = 0; i < tokens.size(); ++i)
      std::cout << i << " : " << tokens[i] << '\n';

#ifdef USE_PROCMON
    FAIL_IF_NOT(tokens[4] == std::string("VmPeak"));
    FAIL_IF_NOT(tokens[9] == std::string("VmRss"));
    FAIL_IF_NOT(tokens[12] == std::string("kB"));
#endif
    // Save the VmRSS size
    vmrss_orig = atof(tokens[10].c_str());
  }

  // Allocate some memory, so we can show that VmRSS increases:
  size_t const N(100000000);
  auto *vdbl = new double[N];
  for (size_t i = 0; i < N; ++i)
    vdbl[i] = std::numbers::pi;

  {
    std::ostringstream msg;

    procmon_resource_print("tst_procmon_use_02", mynode, msg);
    std::string const line = msg.str();
    std::vector<std::string> tokens = rtt_dsxx::tokenize(line, " \t");

    // Save the VmRSS size
    vmrss_bigger = atof(tokens[10].c_str());
  }

  // We should be using more memory.
  FAIL_IF(vmrss_bigger < vmrss_orig);

  for (size_t i = 0; i < N; ++i)
    vdbl[i] = 42.5150;

  delete[] vdbl;

  // Check that the memory use has decreased:
  {
    std::ostringstream msg;
    procmon_resource_print("tst_procmon_use_03", mynode, msg);
    std::string const line = msg.str();
    std::vector<std::string> tokens = rtt_dsxx::tokenize(line, " \t");

    // Save the VmRSS size
    vmrss_end = atof(tokens[10].c_str());
  }

  FAIL_IF(vmrss_bigger < vmrss_end);

  if (ut.numFails == 0)
    PASSMSG("No failures in tst_procmon_basic.");
  return;
}

//------------------------------------------------------------------------------------------------//
void tst_procmon_macro(rtt_dsxx::UnitTest &ut) {
  std::cout << "\nRunning tst_procmon_macro(ut)...\n";

  // if USE_PROCMON is OFF, this will do nothing.  I'm not sure how to test this.
  PROCMON_REPORT("tst_procmon_use_02");

  if (ut.numFails == 0)
    PASSMSG("No failures in tst_procmon_macro.");
  return;
}

} // namespace

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  rtt_dsxx::ScalarUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    tst_procmon_basic(ut);
    tst_procmon_macro(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstProcmon.cc
//------------------------------------------------------------------------------------------------//
