//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   c4/test/tstSwap.cc
 * \author Kent Budge
 * \date   Wed Apr 28 09:31:51 2010
 * \brief  Test c4::determinate_swap and c4::indeterminate_swap functions
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "c4/ParallelUnitTest.hh"
#include "c4/swap.hh"
#include "ds++/Release.hh"
#include <cmath>

using namespace std;
using namespace rtt_dsxx;
using namespace rtt_c4;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

template <typename T> void tstDeterminateSwap(UnitTest &ut) {
  vector<unsigned> outgoing_pid, incoming_pid;
  vector<vector<T>> outgoing_data, incoming_data;

  if (rtt_c4::nodes() == 1) {
    // Should be a no-op.
  } else if (rtt_c4::nodes() == 2) {
    outgoing_pid.resize(1);
    incoming_pid.resize(1);
    if (rtt_c4::node() == 0) {
      outgoing_pid[0] = 1;
      incoming_pid[0] = 1;
    } else {
      Check(rtt_c4::node() == 1);
      outgoing_pid[0] = 0;
      incoming_pid[0] = 0;
    }
  } else {
    Check(rtt_c4::nodes() == 4);
    if (rtt_c4::node() == 0) {
      outgoing_pid.resize(1);
      incoming_pid.resize(1);
      outgoing_pid[0] = 1;
      incoming_pid[0] = 1;
    } else if (rtt_c4::node() == 1) {
      outgoing_pid.resize(2);
      incoming_pid.resize(2);
      outgoing_pid[0] = 0;
      outgoing_pid[1] = 3;
      incoming_pid[0] = 0;
      incoming_pid[1] = 2;
    } else if (rtt_c4::node() == 2) {
      outgoing_pid.resize(1);
      incoming_pid.resize(1);
      outgoing_pid[0] = 1;
      incoming_pid[0] = 3;
    } else {
      outgoing_pid.resize(1);
      incoming_pid.resize(1);
      outgoing_pid[0] = 2;
      incoming_pid[0] = 1;
    }
  }

  outgoing_data.resize(outgoing_pid.size());
  incoming_data.resize(incoming_pid.size());
  for (unsigned i = 0; i < outgoing_pid.size(); ++i) {
    unsigned const pid = outgoing_pid[i];
    outgoing_data[i].resize(2);
    outgoing_data[i][0] = rtt_c4::node();
    outgoing_data[i][1] = pid;
    incoming_data[i].resize(2);
  }
  for (unsigned i = 0; i < incoming_pid.size(); ++i) {
    incoming_data[i].resize(2);
  }
  rtt_c4::determinate_swap(outgoing_pid, outgoing_data, incoming_pid, incoming_data);

  FAIL_IF(incoming_data.size() != incoming_pid.size());
  for (unsigned i = 0; i < incoming_pid.size(); ++i) {
    unsigned const pid = incoming_pid[i];
    FAIL_IF(incoming_data[i].size() != 2);
    FAIL_IF(static_cast<unsigned>(incoming_data[i][0]) != pid ||
            static_cast<unsigned>(incoming_data[i][1]) != static_cast<unsigned>(rtt_c4::node()));
  }

  // Second version (no processor list)

  outgoing_data.resize(0);
  outgoing_data.resize(rtt_c4::nodes());
  incoming_data.resize(0);
  incoming_data.resize(rtt_c4::nodes());
  for (unsigned int pid : outgoing_pid) {
    outgoing_data[pid].resize(2);
    outgoing_data[pid][0] = rtt_c4::node();
    outgoing_data[pid][1] = pid;
  }
  for (unsigned int pid : incoming_pid) {
    incoming_data[pid].resize(2);
  }
  rtt_c4::determinate_swap(outgoing_data, incoming_data);

  for (unsigned int pid : incoming_pid) {
    FAIL_IF(incoming_data[pid].size() != 2);
    FAIL_IF(static_cast<unsigned>(incoming_data[pid][0]) != pid ||
            static_cast<unsigned>(incoming_data[pid][1]) != static_cast<unsigned>(rtt_c4::node()));
  }
}

void tstSemideterminateSwap(UnitTest &ut) {
  vector<unsigned> outgoing_pid, incoming_pid;
  vector<vector<unsigned>> outgoing_data, incoming_data;

  if (rtt_c4::nodes() == 1) {
    // Should be a no-op.
  } else if (rtt_c4::nodes() == 2) {
    outgoing_pid.resize(1);
    incoming_pid.resize(1);
    if (rtt_c4::node() == 0) {
      outgoing_pid[0] = 1;
      incoming_pid[0] = 1;
    } else {
      Check(rtt_c4::node() == 1);
      outgoing_pid[0] = 0;
      incoming_pid[0] = 0;
    }
  } else {
    Check(rtt_c4::nodes() == 4);
    if (rtt_c4::node() == 0) {
      outgoing_pid.resize(1);
      incoming_pid.resize(1);
      outgoing_pid[0] = 1;
      incoming_pid[0] = 1;
    } else if (rtt_c4::node() == 1) {
      outgoing_pid.resize(2);
      incoming_pid.resize(2);
      outgoing_pid[0] = 0;
      outgoing_pid[1] = 3;
      incoming_pid[0] = 0;
      incoming_pid[1] = 2;
    } else if (rtt_c4::node() == 2) {
      outgoing_pid.resize(1);
      incoming_pid.resize(1);
      outgoing_pid[0] = 1;
      incoming_pid[0] = 3;
    } else {
      outgoing_pid.resize(1);
      incoming_pid.resize(1);
      outgoing_pid[0] = 2;
      incoming_pid[0] = 1;
    }
  }

  outgoing_data.resize(outgoing_pid.size());
  incoming_data.resize(incoming_pid.size());
  for (unsigned i = 0; i < outgoing_pid.size(); ++i) {
    unsigned const pid = outgoing_pid[i];
    outgoing_data[i].resize(2);
    outgoing_data[i][0] = rtt_c4::node();
    outgoing_data[i][1] = pid;
    //        incoming_data[i].resize(2);
  }
  for (unsigned i = 0; i < incoming_pid.size(); ++i) {
    incoming_data[i].resize(2);
  }
  rtt_c4::semideterminate_swap(outgoing_pid, outgoing_data, incoming_pid, incoming_data);

  if (incoming_data.size() != incoming_pid.size()) {
    ut.failure("Incoming data is NOT correct count");
  } else {
    ut.passes("Incoming data is correct count");
  }
  for (unsigned i = 0; i < incoming_pid.size(); ++i) {
    unsigned const pid = incoming_pid[i];
    if (incoming_data[i].size() != 2) {
      ut.failure("Incoming data is NOT correct size");
    } else {
      ut.passes("Incoming data is correct size");
    }
    if (incoming_data[i][0] != pid || incoming_data[i][1] != static_cast<size_t>(rtt_c4::node())) {
      ut.failure("Incoming data is NOT correct");
    } else {
      ut.passes("Incoming data is correct");
    }
  }
}

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  rtt_c4::ParallelUnitTest ut(argc, argv, release);
  try {
    tstDeterminateSwap<unsigned>(ut);
    tstDeterminateSwap<double>(ut);
    tstSemideterminateSwap(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstSwap.cc
//------------------------------------------------------------------------------------------------//
