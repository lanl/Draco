//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   c4/test/tstopstream.cc
 * \author Kent Budge
 * \date   Wed Apr 28 09:31:51 2010
 * \brief  Test c4::determinate_swap and c4::indeterminate_swap functions
 * \note   Copyright (C) 2018-2022 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

#include "c4/ParallelUnitTest.hh"
#include "c4/opstream.hh"
#include "ds++/Release.hh"
#include <cmath>

using namespace std;
using namespace rtt_dsxx;
using namespace rtt_c4;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

void tstopstream(UnitTest &ut) {
  unsigned const pid = rtt_c4::node();
  if (pid == 0)
    cout << "Start of write:" << endl;

  opstream pout;

  pout << "MPI rank " << pid << " reporting ..." << endl;

  pout.send();

  pout.shrink_to_fit();

  pout << "MPI rank " << pid << " reporting a second time ..." << endl;

  pout.shrink_to_fit();
  pout.send();

  if (pid == 0)
    cout << ": End of write" << endl;

  {
    // Test dynamic object creation and destruction
    auto *ptr_pout = new opstream;
    delete ptr_pout;
  }

  PASSMSG("completed serialized write without hanging or segfaulting");
}

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
