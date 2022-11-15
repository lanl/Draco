//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   c4/test/tstBroadcast.cc
 * \author Thomas M. Evans
 * \date   Tue Apr  2 15:57:11 2002
 * \brief  c4 Broadcast communication test.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "c4/ParallelUnitTest.hh"
#include "ds++/Release.hh"
#include "ds++/Soft_Equivalence.hh"
#include <sstream>

using namespace std;

using rtt_c4::broadcast;
using rtt_dsxx::soft_equiv;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

void test_simple(rtt_dsxx::UnitTest &ut) {

  char c = 0;
  unsigned char uc = 0;
  int i = 0;
  long l = 0;
  long long ll = 0;
  unsigned long long ull = 0;
  unsigned long ul = 0;
  unsigned short us = 0;
  short s = 0;
  float f = 0;
  double d = 0;
  long double ld = 0;
  vector<double> vref(10, 3.1415);
  vector<double> v(10, 0.0);
  vector<float> vempty_ref;
  vector<float> vempty;
  string msgref("hello, world!");
  string msg;

  // assign on node 0
  if (rtt_c4::node() == 0) {
    c = 'A';
    uc = 'B';
    i = 1;
    l = 1000;
    ll = 1000;
    ull = 1000;
    ul = 1000;
    us = 1000;
    s = 1000;
    f = 1.5F;
    d = 2.5;
    ld = 3.5L;
    v = vref;
    msg = msgref;
  } else {
    // reserve enough space to receive the broadcast string.
    msg.resize(msgref.length());
  }

  // send out data, using node 0 as root
  broadcast(&c, 1, 0);
  broadcast(&uc, 1, 0);
  broadcast(&i, 1, 0);
  broadcast(&l, 1, 0);
  broadcast(&ll, 1, 0);
  broadcast(&ull, 1, 0);
  broadcast(&ul, 1, 0);
  broadcast(&us, 1, 0);
  broadcast(&s, 1, 0);
  broadcast(&f, 1, 0);
  broadcast(&d, 1, 0);
#ifndef CRAYPE_CCE
  broadcast(&ld, 1, 0);
#endif
  broadcast(v.begin(), v.end(), v.begin());
  broadcast(vempty.begin(), vempty.end(), vempty.begin());
  broadcast(msg.begin(), msg.end(), msg.begin());

  // check scalar values
  FAIL_IF_NOT(c == 'A');
  FAIL_IF_NOT(uc == 'B');
  FAIL_IF_NOT(i == 1);
  FAIL_IF_NOT(l == 1000);
  FAIL_IF_NOT(ll == 1000);
  FAIL_IF_NOT(ull == 1000);
  FAIL_IF_NOT(ul == 1000);
  FAIL_IF_NOT(us == 1000);
  FAIL_IF_NOT(s == 1000);
  FAIL_IF_NOT(soft_equiv(f, 1.5F));
  FAIL_IF_NOT(soft_equiv(d, 2.5));
#ifndef CRAYPE_CCE
  FAIL_IF_NOT(soft_equiv(ld, 3.5L));
#endif
  // check array values
  FAIL_IF_NOT(soft_equiv(v.begin(), v.end(), vref.begin(), vref.end()));
  FAIL_IF_NOT(
      soft_equiv(vempty.begin(), vempty.end(), vempty_ref.begin(), vempty_ref.end(), 1.0e-7F));
  FAIL_IF_NOT(msg == msgref);

  // safer 4 argument form (in case msg has not been resized).
  if (rtt_c4::node() != 0)
    msg = "foo bar baz 9"; // same length as msgref.
  broadcast(msg.begin(), msg.end(), msg.begin(), msg.end());
  FAIL_IF_NOT(msg == msgref);

  try {
    string badmsg; // length never set.
    broadcast(msg.begin(), msg.end(), badmsg.begin(), badmsg.end());
    // The above command should throw on all processors.
    FAIL_IF_NOT(rtt_c4::node() == 0);
  } catch (exception & /*error*/) {
    ostringstream mymsg;
    mymsg << "Successfully caught a range violation in broadcast on PE " << rtt_c4::node();
    PASSMSG(mymsg.str());
  }

  rtt_c4::global_barrier();

  if (ut.numFails == 0) {
    ostringstream m;
    m << "test_simple() okay on " << rtt_c4::node();
    PASSMSG(m.str());
  }
  return;
}

//------------------------------------------------------------------------------------------------//
//  By adjusting the parameters below, this test will overflow the MPI memory
//  buffers.  Read the comments below if you'd like to do this.
void test_loop(rtt_dsxx::UnitTest &ut) {
  // save state
  unsigned const nf(ut.numFails);

  // kmax controls how much data is broadcast.  If kmax is too big (like
  // 10000000), shmem will fail.
  int const kmax = 10;

  if (rtt_c4::node() == 0) // host proc
  {
    // send out the values on host
    for (int k = 0; k < kmax; ++k) {
      Insist(!broadcast(&k, 1, 0), "MPI Error");
      double foo = k + 0.5;
      Insist(!broadcast(&foo, 1, 0), "MPI Error");
    }
  } else // all other processors
  {
    // Use sleep() if you want the host processor to fill up the buffers.  We
    // comment out the sleep() command here because it's not supported on all
    // platforms.

    // sleep(10);

    int kk = 0;
    double foofoo = 0.0;
    for (int k = 0; k < kmax; ++k) {
      kk = -1;
      foofoo = -2.0;
      Insist(!broadcast(&kk, 1, 0), "MPI Error");
      FAIL_IF_NOT(kk == k);
      Insist(!broadcast(&foofoo, 1, 0), "MPI Error");
      FAIL_IF_NOT(rtt_dsxx::soft_equiv(foofoo, k + 0.5));
    }
  }

  if (ut.numFails == nf) {
    ostringstream m;
    m << "test_loop() okay on " << rtt_c4::node();
    PASSMSG(m.str());
  }
  return;
}

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  rtt_c4::ParallelUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    // Unit tests
    test_simple(ut);
    test_loop(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstBroadcast.cc
//------------------------------------------------------------------------------------------------//
