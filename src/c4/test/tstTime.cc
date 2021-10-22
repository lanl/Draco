//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   c4/test/tstTime.cc
 * \author Thomas M. Evans
 * \date   Mon Mar 25 17:19:16 2002
 * \brief  Test timing functions in C4.
 * \note   Copyright (C) 2010-2021 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "c4/Global_Timer.hh"
#include "c4/ParallelUnitTest.hh"
#include "ds++/Release.hh"
#include "ds++/Soft_Equivalence.hh"
#include <sstream>

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

void wall_clock_test(rtt_dsxx::UnitTest &ut) {
  using std::cout;
  using std::endl;
  using std::ostringstream;
  using std::set;
  using std::string;

  using rtt_c4::Global_Timer;
  using rtt_c4::Timer;
  using rtt_c4::wall_clock_resolution;
  using rtt_c4::wall_clock_time;
  using rtt_dsxx::soft_equiv;

  Global_Timer do_timer("do_timer");
  Global_Timer do_not_timer("do_not_timer");

  if (do_timer.is_active()) {
    ut.failure("do_timer should NOT be active");
  }

  set<string> active_timers;
  active_timers.insert("do_timer");
  active_timers.insert("do_global_timer");
  Global_Timer::set_selected_activity(active_timers, true);

  if (rtt_c4::nodes() == 1 && !do_timer.is_active()) {
    ut.failure("do_timer SHOULD be active");
  }

  do_timer.start();

  double const wcr(rtt_c4::wall_clock_resolution());
  if (wcr > 0.0 && wcr <= 1000.0) {
    ostringstream msg;
    msg << "The timer has a wall clock resoution of " << wcr << " ticks." << endl;
    PASSMSG(msg.str());
  } else {
    ostringstream msg;
    msg << "The timer does not appear to have a reasonable resolution."
        << " rtt_c4::wall_clock_resolution() = " << wcr << " ticks." << endl;
    FAILMSG(msg.str());
  }

  Timer t;

#ifdef _MSC_VER
  // KT: I'm not sure why I'm having so much trouble on Win32 for the timing precision check.  I may
  // need to use a different mechanism for timing values.  Right now the ngihtly tests fail about 1
  // out of 20 times.  When the check fails it is significant:
  //
  // t.wall_clock() value does not match the expected value.
  //   end            = 43461.5
  //   begin          = 43460.1
  //   end-begin      = 1.39073
  //   t.wall_clock() = 1.40682
  //   error          = 0.0160977
  //   prec           = 0.002
  //
  // This may be due to the system being busy or scans done by LANL IT.  I'm not really sure what is
  // going on.  It always passes for interactive jobs.
  //
  // In any case, I am making the prec value for Win32 10x the Linux value.  This should keep the
  // tests passing and is still valid because the relative error is less than 10% for this very
  // short time interval.
  double const prec(20.0 * t.posix_err());
#else
  double const prec(4.0 * t.posix_err());
#endif
  double begin(rtt_c4::wall_clock_time());

  t.start();

  FAIL_IF_NOT(t.on());

  // do some work
  if (rtt_c4::node() == 0)
    std::cout << "\nDoing some work..." << std::endl;
  size_t len(10000000);
  std::vector<double> foo(len);
  double sum(0);
  for (size_t i = 0; i < len; ++i) {
    auto const d(static_cast<double>(i + 1));
    foo[i] = std::sqrt(std::log(d * 3.14) * std::fabs(std::cos(d / 3.14)));
    sum += foo[i];
  }

  double end = rtt_c4::wall_clock_time();
  t.stop();

  FAIL_IF(t.on());

  double const error(t.wall_clock() - (end - begin));
  if (std::fabs(error) <= prec) {
    PASSMSG("wall_clock() value looks ok.");
  } else {
    ostringstream msg;
    msg << "t.wall_clock() value does not match the expected value."
        << "\n\tend            = " << end << "\n\tbegin          = " << begin
        << "\n\tend-begin      = " << end - begin << "\n\tt.wall_clock() = " << t.wall_clock()
        << "\n\terror          = " << error << "\n\tprec           = " << prec << endl;
    FAILMSG(msg.str());
  }

  //----------------------------------------------------------------------------------------------//
  // Ensure that system + user <= wall
  //
  // Due to round off errors, the wall clock time might be less than the system + user time.  But
  // this difference should never exceed t.posix_err().
  //----------------------------------------------------------------------------------------------//

  double const deltaWallTime(t.wall_clock() - (t.system_cpu() + t.user_cpu()));
#ifdef _MSC_VER
  double const time_resolution(1.0);
#else
  double const time_resolution(prec);
#endif
  if (deltaWallTime > 0.0 || std::fabs(deltaWallTime) <= time_resolution) {
    ostringstream msg;
    msg << "The sum of cpu and user time is less than or equal to the\n\t"
        << "reported wall clock time (within error bars = " << time_resolution << " secs.)."
        << endl;
    PASSMSG(msg.str());
  } else {
    ostringstream msg;
    msg << "The sum of cpu and user time exceeds the reported wall "
        << "clock time.  Here are the details:"
        << "\n\tposix_error() = " << prec << " sec."
        << "\n\tdeltaWallTime = " << deltaWallTime << " sec."
        << "\n\tSystem time   = " << t.system_cpu() << " sec."
        << "\n\tUser time     = " << t.user_cpu() << " sec."
        << "\n\tWall time     = " << t.wall_clock() << " sec." << endl;
    FAILMSG(msg.str());
  }

  //------------------------------------------------------//
  // Demonstrate print functions:
  //------------------------------------------------------//

  if (rtt_c4::node() == 0)
    cout << "Demonstration of the print() member function via the\n"
         << "\toperator<<(ostream&,Timer&) overloaded operator.\n"
         << endl;

  cout << "Timer = " << t << endl;

  //------------------------------------------------------//
  // Do a second timing:
  //------------------------------------------------------//

  cout << "\nCreate a Timer Report after two timing cycles:\n" << endl;

  t.start();
  for (size_t i = 0; i < len; ++i)
    foo[i] = static_cast<double>(i * 4);
  t.stop();

  t.print(cout, 6);

  // Test the single line printout
  if (rtt_c4::node() == 0) {
    std::ostringstream timingsingleline;
    t.printline(timingsingleline, 4, 8);
#ifdef HAVE_PAPI
    if (timingsingleline.str().length() == 34)
#else
    if (timingsingleline.str().length() == 26)
#endif
      PASSMSG(string("printline() returned a single line of the") + " expected length.");
    else
      FAILMSG(string("printline() did not return a line of the ") + "expected length.");
  }

  //------------------------------------------------------//
  // Check the number of intervals
  //------------------------------------------------------//

  int const expectedNumberOfIntervals(2);
  if (t.intervals() == expectedNumberOfIntervals)
    PASSMSG("Found the expected number of intervals.");
  else
    FAILMSG("Did not find the expected number of intervals.");

  //------------------------------------------------------//
  // Check the merge method
  //------------------------------------------------------//

  double const old_wall_time = t.sum_wall_clock();
  double const old_system_time = t.sum_system_cpu();
  double const old_user_time = t.sum_user_cpu();
  int const old_intervals = t.intervals();
  t.merge(t);

  if (rtt_dsxx::soft_equiv(2 * old_wall_time, t.sum_wall_clock()) &&
      rtt_dsxx::soft_equiv(2 * old_system_time, t.sum_system_cpu()) &&
      rtt_dsxx::soft_equiv(2 * old_user_time, t.sum_user_cpu()) &&
      2 * old_intervals == t.intervals())
    PASSMSG("merge okay");
  else
    FAILMSG("merge NOT okay");

  //------------------------------------------------------------//
  // Check PAPI data
  //------------------------------------------------------------//

  long long cachemisses = t.sum_cache_misses();
  long long cachehits = t.sum_cache_hits();
  long long flops = t.sum_floating_operations();

#ifdef HAVE_PAPI

  std::cout << "PAPI metrics report:\n"
            << "   Cache misses : " << cachemisses << "\n"
            << "   Cache hits   : " << cachehits << "\n"
            << "   FLOP         : " << flops << std::endl;

  if (cachemisses == 0 && cachehits == 0 && flops == 0)
    FAILMSG("PAPI metrics returned 0 when PAPI was available.");
  else
    PASSMSG("PAPI metrics returned >0 values when PAPI is available.");

#else
  if (cachemisses == 0 && cachehits == 0 && flops == 0)
    PASSMSG("PAPI metrics return 0 when PAPI is not available.");
  else
    FAILMSG("PAPI metrics did not return 0 when PAPI was not available.");
#endif

  // Exercise the global report and reset. Unfortunately, there's no easy way to
  // check the validity of the output except to eyeball.

  do_timer.stop();
  Global_Timer::report_all(cout);
  Global_Timer::reset_all();
  Global_Timer::report_all(cout);

  return;
}

//------------------------------------------------------------------------------------------------//
void test_pause(rtt_dsxx::UnitTest &ut) {
  using rtt_c4::Timer;

  std::cout << "Starting tstTime::test_pause tests..." << std::endl;

  double begin(rtt_c4::wall_clock_time());

  double const pauseTime(1.0); // seconds

  Timer::pause(pauseTime);

  double end(rtt_c4::wall_clock_time());

  FAIL_IF(end - begin < pauseTime);

  std::cout << "Starting tstTime::test_pause tests...done\n" << std::endl;

  return;
}

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  rtt_c4::ParallelUnitTest ut(argc, argv, rtt_dsxx::release);
  rtt_c4::Timer::initialize(argc, argv);
  try {
    wall_clock_test(ut);
    test_pause(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstTime.cc
//------------------------------------------------------------------------------------------------//
