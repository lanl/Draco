//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   diagnostics/test/tstTiming.cc
 * \author Thomas M. Evans
 * \date   Mon Dec 12 15:32:10 2005
 * \brief  Test the diagnostics/TIMER macros
 * \note   Copyright (C) 2010-2023 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "c4/ParallelUnitTest.hh"
#include "diagnostics/Diagnostics.hh"
#include "diagnostics/Timing.hh"
#include "ds++/Release.hh"
#include "ds++/Soft_Equivalence.hh"
#include <iomanip>

using namespace std;
using rtt_dsxx::soft_equiv;
using D = rtt_diagnostics::Timing_Diagnostics;

//------------------------------------------------------------------------------------------------//
// TEST HELPERS
//------------------------------------------------------------------------------------------------//

void do_A() {
  TIMER(A_timer);
  TIMER_START("A_iteration", A_timer);

  // do a mat-vec multiply

  int S = 9000;

  vector<double> b(S, 0.0);
  vector<double> x(S, 0.0);

  double A = 1.0;
  double B = 1.1;
  double C = 1.01;

  x[0] = 0.2;
  for (int i = 1; i < S; i++)
    x[i] = i + 1.1 + x[i - 1];

  for (int i = 1; i < (S - 1); ++i) {
    b[i] = x[i - 1] * A + x[i] * B + x[i + 1] * C;
  }

  b[0] = B * x[0] + C * x[1];
  b[S - 1] = A * x[S - 2] + B * x[S - 1];

  TIMER_STOP("A_iteration", A_timer);
  TIMER_RECORD("A_iteration", A_timer);
}

//------------------------------------------------------------------------------------------------//

void do_B() {
  TIMER(B_timer);
  TIMER_START("B_iteration", B_timer);

  // do a mat-vec multiply

  int S = 6000;

  vector<double> b(S, 0.0);
  vector<double> x(S, 0.0);

  double A = 1.0;
  double B = 1.1;
  double C = 1.01;

  x[0] = 0.2;
  for (int i = 1; i < S; i++)
    x[i] = i + 1.1 + x[i - 1];

  for (int i = 1; i < (S - 1); ++i) {
    b[i] = x[i - 1] * A + x[i] * B + x[i + 1] * C;
  }

  b[0] = B * x[0] + C * x[1];
  b[S - 1] = A * x[S - 2] + B * x[S - 1];

  TIMER_STOP("B_iteration", B_timer);
  TIMER_RECORD("B_iteration", B_timer);
}

//------------------------------------------------------------------------------------------------//

void do_C() {
  TIMER(C_timer);
  TIMER_START("C_iteration", C_timer);

  // do a mat-vec multiply

  int S = 3000;

  vector<double> b(S, 0.0);
  vector<double> x(S, 0.0);

  double A = 1.0;
  double B = 1.1;
  double C = 1.01;

  x[0] = 0.2;
  for (int i = 1; i < S; i++)
    x[i] = i + 1.1 + x[i - 1];

  for (int i = 1; i < (S - 1); ++i) {
    b[i] = x[i - 1] * A + x[i] * B + x[i + 1] * C;
  }

  b[0] = B * x[0] + C * x[1];
  b[S - 1] = A * x[S - 2] + B * x[S - 1];

  TIMER_STOP("C_iteration", C_timer);
  TIMER_RECORD("C_iteration", C_timer);
}

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

void timing_active() {
#ifdef DRACO_TIMING_ON
  cout << ">>> Testing timing macros with value " << DRACO_TIMING << endl;
#else
  cout << ">>> Timing macros inactive" << endl;
#endif
}

//------------------------------------------------------------------------------------------------//

void test_timing(rtt_dsxx::UnitTest &ut) {
  // add to some timers
  D::update_timer("A", 1.2);
  D::update_timer("B", 1.1);
  D::update_timer("B", 2.3);

  FAIL_IF_NOT(soft_equiv(D::timer_value("A"), 1.2));
  FAIL_IF_NOT(soft_equiv(D::timer_value("B"), 3.4));
  FAIL_IF_NOT(soft_equiv(D::timer_value("C"), 0.0));

  D::reset_timer("B");
  D::update_timer("A", 1.3);
  FAIL_IF_NOT(soft_equiv(D::timer_value("A"), 2.5));
  FAIL_IF_NOT(soft_equiv(D::timer_value("B"), 0.0));
  FAIL_IF_NOT(soft_equiv(D::timer_value("C"), 0.0));

  vector<string> timers = D::timer_keys();
  FAIL_IF_NOT(timers.size() == 3);
  FAIL_IF_NOT(timers[0] == "A");
  FAIL_IF_NOT(timers[1] == "B");
  FAIL_IF_NOT(timers[2] == "C");

  D::delete_timer("B");
  timers = D::timer_keys();
  FAIL_IF_NOT(timers.size() == 2);
  FAIL_IF_NOT(timers[0] == "A");
  FAIL_IF_NOT(timers[1] == "C");

  // calling timer_value on B will get it back
  FAIL_IF_NOT(soft_equiv(D::timer_value("A"), 2.5));
  FAIL_IF_NOT(soft_equiv(D::timer_value("B"), 0.0));
  FAIL_IF_NOT(soft_equiv(D::timer_value("C"), 0.0));
  timers = D::timer_keys();
  FAIL_IF_NOT(timers.size() == 3);
  FAIL_IF_NOT(D::num_timers() == 3);

  // delete all timers
  D::delete_timers();
  FAIL_IF_NOT(D::num_timers() == 0);
  timers = D::timer_keys();
  FAIL_IF_NOT(timers.size() == 0);

  D::update_timer("B", 12.4);
  D::update_timer("C", 1.3);
  FAIL_IF_NOT(soft_equiv(D::timer_value("A"), 0.0));
  FAIL_IF_NOT(soft_equiv(D::timer_value("B"), 12.4));
  FAIL_IF_NOT(soft_equiv(D::timer_value("C"), 1.3));

  // reset all timers
  D::reset_timers();
  FAIL_IF_NOT(soft_equiv(D::timer_value("A"), 0.0));
  FAIL_IF_NOT(soft_equiv(D::timer_value("B"), 0.0));
  FAIL_IF_NOT(soft_equiv(D::timer_value("C"), 0.0));

  if (ut.numFails == 0)
    PASSMSG("Diagnostics timer lists ok.");
}

//------------------------------------------------------------------------------------------------//
void test_macros(rtt_dsxx::UnitTest &ut) {
  D::delete_timers(); // delete all existing timers

  // make timers and do results
  TIMER(outer_timer);
  TIMER_START("Outer", outer_timer);
  do_A();
  do_B();
  do_C();

  TIMER_STOP("Outer", outer_timer);
  TIMER_RECORD("Outer", outer_timer);

  // if the timers are off we get no timing data
  vector<string> keys = D::timer_keys();
#if DRACO_TIMING == 0
  FAIL_IF_NOT(keys.size() == 0);
  FAIL_IF_NOT(D::num_timers() == 0);
#else
#ifndef DRACO_CALIPER
  FAIL_IF(keys.size() != 4);
  cout << setw(15) << "Routine" << setw(15) << "Fraction\n" << string(30, '-') << endl;

  // get the keys and print a table
  double const total = D::timer_value("Outer");
  FAIL_IF(total <= 0.0);

  cout.precision(4);
  cout.setf(ios::fixed, ios::floatfield);

  for (auto const &key : keys) {
    double const fraction = D::timer_value(key) / total;
    cout << setw(15) << key << setw(15) << fraction << endl;
  }
  cout << "The total time was " << total << "\n" << endl;
#endif // DRACO_CALIPER
#endif

  TIMER_REPORT(outer_timer, cout, "Total time");

  if (ut.numFails == 0)
    PASSMSG("Timer macros ok.");
}

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  rtt_c4::ParallelUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    timing_active();
    test_timing(ut);
    test_macros(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstTiming.cc
//------------------------------------------------------------------------------------------------//
