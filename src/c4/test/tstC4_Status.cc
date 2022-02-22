//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   c4/test/tstC4_Status.cc
 * \author Robert B. Lowrie
 * \date   Friday May 26 19:58:19 2017
 * \brief  Unit test for C4_Status class.
 * \note   Copyright (C) 2017-2022 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

#include "c4/ParallelUnitTest.hh"
#include "ds++/Release.hh"
#include <vector>

using namespace std;
using rtt_c4::C4_Req;
using rtt_c4::C4_Status;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

void tst2Procs(rtt_dsxx::UnitTest &ut) {

  PASSMSG("Running tst2Procs.");

  const int my_proc = rtt_c4::node();

  {
    const uint32_t my_nranks = rtt_c4::nranks();
    const int my_nodes = rtt_c4::nodes();
    FAIL_IF_NOT(static_cast<uint32_t>(my_nodes) == my_nranks);
    FAIL_IF_NOT(rtt_c4::isMpiInit());
  }

  C4_Status status{};
  C4_Req request;

  const int num_int = 2;
  const int num_double = 5;
  const int tag = 101;

  if (my_proc == 0) {
    vector<int> send_buffer(num_int);
    vector<double> recv_buffer(num_double);
    rtt_c4::receive_async(request, &recv_buffer[0], num_double, rtt_c4::any_source, tag);
    rtt_c4::send_async(&send_buffer[0], num_int, 1, tag);
    request.wait(&status);
    if (status.get_source() == 1)
      PASSMSG("get_source() passed on processor 0");
    else
      FAILMSG("get_source() failed on processor 0");
    if (status.get_message_size() == num_double * sizeof(double))
      PASSMSG("get_message_size() passed on processor 0");
    else
      FAILMSG("get_message_size() failed on processor 0");
    if (status.get_status_obj())
      PASSMSG("get_status_obj() passed on processor 0");
    else
      FAILMSG("get_status_obj() failed on processor 0");
  } else { // my_proc == 1
    vector<double> send_buffer(num_double);
    vector<int> recv_buffer(num_int);
    rtt_c4::receive_async(request, &recv_buffer[0], num_int, rtt_c4::any_source, tag);
    rtt_c4::send_async(&send_buffer[0], num_double, 0, tag);
    request.wait(&status);
    if (status.get_source() == 0)
      PASSMSG("get_source() passed on processor 1");
    else
      FAILMSG("get_source() failed on processor 1");
    if (status.get_message_size() == num_int * sizeof(int))
      PASSMSG("get_message_size() passed on processor 1");
    else
      FAILMSG("get_message_size() failed on processor 1");
    if (status.get_status_obj())
      PASSMSG("get_status_obj() passed on processor 1");
    else
      FAILMSG("get_status_obj() failed on processor 1");
  }

  if (ut.numFails == 0)
    PASSMSG("tstC4_Status() is okay.");

  return;
}

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  rtt_c4::ParallelUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    if (rtt_c4::nodes() == 2)
      tst2Procs(ut);
    else
      FAILMSG("tstC4_Status should only be run on 2 processors!");
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstC4_Status.cc
//------------------------------------------------------------------------------------------------//
