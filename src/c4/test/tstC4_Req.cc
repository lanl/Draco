//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   c4/test/tstC4_Req.cc
 * \author Kelly Thompson
 * \date   Tue Nov  1 15:49:44 2005
 * \brief  Unit test for C4_Req class.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "c4/ParallelUnitTest.hh"
#include "ds++/Release.hh"
#include <array>

using namespace std;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

void tstCopyConstructor(rtt_dsxx::UnitTest &ut) {
  using rtt_c4::C4_Req;

  C4_Req requestA;
  C4_Req requestB(requestA);

  // The behavior of the copy constructor is not obvious.  If requestA has not
  // been used (inuse() returns false) then requestA != requestB.

  if (!requestA.inuse() && requestA == requestB)
    FAILMSG("requestA.inuse() is false, so requestA cannot == requestB.");

  if (!requestA.inuse() && requestA != requestB)
    PASSMSG("requestA.inuse() is false and requestA != requestB.");

  if (requestA.inuse() && requestA == requestB)
    PASSMSG("requestA.inuse() is true and requestA == requestB.");

  if (requestA.inuse() && requestA != requestB)
    FAILMSG("requestA.inuse() is true, so requestA must == requestB.");

  if (ut.numFails == 0)
    PASSMSG("tstCopyConstructor() is okay.");

  return;
}

//------------------------------------------------------------------------------------------------//
void tstTraits(rtt_dsxx::UnitTest &ut) {
  using rtt_c4::C4_Traits;

  {
    constexpr bool is_uchar = C4_Traits<unsigned char>::tag == 432;
    FAIL_IF_NOT(is_uchar);
    constexpr bool is_short = C4_Traits<short>::tag == 433;
    FAIL_IF_NOT(is_short);
    constexpr bool is_ushort = C4_Traits<unsigned short>::tag == 434;
    FAIL_IF_NOT(is_ushort);
    constexpr bool is_uint = C4_Traits<unsigned int>::tag == 436;
    FAIL_IF_NOT(is_uint);
    constexpr bool is_ulong = C4_Traits<unsigned long>::tag == 438;
    FAIL_IF_NOT(is_ulong);
    constexpr bool is_longdouble = C4_Traits<long double>::tag == 441;
    FAIL_IF_NOT(is_longdouble);
  }
#ifdef C4_MPI
  {
    using rtt_c4::MPI_Traits;
    FAIL_IF_NOT(MPI_Traits<unsigned char>::element_type() == MPI_UNSIGNED_CHAR);
    FAIL_IF_NOT(MPI_Traits<short>::element_type() == MPI_SHORT);
    FAIL_IF_NOT(MPI_Traits<unsigned short>::element_type() == MPI_UNSIGNED_SHORT);
    FAIL_IF_NOT(MPI_Traits<unsigned int>::element_type() == MPI_UNSIGNED);
    FAIL_IF_NOT(MPI_Traits<unsigned long>::element_type() == MPI_UNSIGNED_LONG);
    FAIL_IF_NOT(MPI_Traits<long double>::element_type() == MPI_LONG_DOUBLE);
  }
#endif

  return;
}

//------------------------------------------------------------------------------------------------//
void tstWait(rtt_dsxx::UnitTest &ut) {
  using namespace rtt_c4;

  if (rtt_c4::node() > 0) {
    cout << "sending from processor " << get_processor_name() << ':' << endl;
    array<int, 1> buffer{0};
    buffer[0] = node();
    C4_Req outgoing = send_async(buffer.data(), 1U, 0);
    unsigned result = wait_any(1U, &outgoing);
    FAIL_IF_NOT(result == 0);
  } else {
    cout << "receiving to processor " << get_processor_name() << ':' << endl;
    Check(rtt_c4::nodes() < 5);
    array<C4_Req, 4> requests;
    array<bool, 4> done = {false, false, false, false};
    for (int p = 1; p < nodes(); ++p) {
      array<array<int, 1>, 4> buffer{};
      requests[p] = receive_async(buffer[p].data(), 1U, p);
      done[p] = false;
    }
    for (int c = 1; c < nodes(); ++c) {
      unsigned result = wait_any(nodes(), requests.data());
      FAIL_IF(done[result]);
      done[result] = true;
    }
  }
  return;
}

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  rtt_c4::ParallelUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    tstCopyConstructor(ut);
    tstTraits(ut);
    tstWait(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstC4_Req.cc
//------------------------------------------------------------------------------------------------//
