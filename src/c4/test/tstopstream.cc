//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   c4/test/tstopstream.cc
 * \author Kent Budge
 * \brief  Test c4::determinate_swap and c4::indeterminate_swap functions
 * \note   Copyright (C) 2018-2026 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "dsxx/Release.hh"
#include "dsxx/UnitTest.hh"
#include "c4/C4_Functions.hh"
#include "c4/ParallelUnitTest.hh"
#include "c4/opstream.hh"
#include <iostream>

namespace {

using namespace std;
using namespace rtt_dsxx;
using namespace rtt_c4;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

void tstopstream(UnitTest &ut) {
  unsigned const pid = rtt_c4::node();
  if (pid == 0)
    cout << "Start of write:\n";

  opstream pout;

  pout << "MPI rank " << pid << " reporting ...\n";

  pout.send();
  pout.shrink_to_fit();
  pout << "MPI rank " << pid << " reporting a second time ...\n";
  pout.shrink_to_fit();
  pout.send();

  if (pid == 0)
    cout << ": End of write\n";

  {
    // Test dynamic object creation and destruction
    auto *ptr_pout = new opstream;
    delete ptr_pout;
  }

  PASSMSG("completed serialized write without hanging or segfaulting");
}

} // namespace

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  rtt_c4::ParallelUnitTest ut(argc, argv, release);
  try {
    tstopstream(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstopstream.cc
//------------------------------------------------------------------------------------------------//
