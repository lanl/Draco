//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   c4/Timer.cc
 * \author Thomas M. Evans
 * \date   Mon Mar 25 17:56:11 2002
 * \note   Copyright (C) 2010-2023 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "Timer.hh"
#include "ds++/XGetopt.hh"
#include <array>
#include <cmath>
#include <cstdlib>
#include <iomanip>

namespace rtt_c4 {

//------------------------------------------------------------------------------------------------//
// Constructor
//------------------------------------------------------------------------------------------------//

//! Constructor
Timer::Timer()
    : tms_begin(DRACO_TIME_TYPE()), tms_end(DRACO_TIME_TYPE()),
      posix_clock_ticks_per_second(static_cast<int>(DRACO_CLOCKS_PER_SEC)),
      isMPIWtimeAvailable(setIsMPIWtimeAvailable()) {
#if defined(_MSC_VER)
  static_assert(DRACO_CLOCKS_PER_SEC < INT32_MAX, "!(DRACO_CLOCKS_PER_SEC < INT32_MAX)");
#else
  Check(DRACO_CLOCKS_PER_SEC < INT32_MAX);
#endif
  reset();
}

//------------------------------------------------------------------------------------------------//
// Member functions
//------------------------------------------------------------------------------------------------//

//! Print out a timing report.
void Timer::print(std::ostream &out, int p) const {
  using std::ios;
  using std::setw;

  out.setf(ios::fixed, ios::floatfield);
  out.precision(p);
  out << '\n';

  if (num_intervals > 1)
    out << "LAST INTERVAL: " << '\n';

  out << setw(20) << "WALL CLOCK TIME: " << wall_clock() << " sec." << '\n'
      << setw(20) << "  USER CPU TIME: " << user_cpu() << " sec." << '\n'
      << setw(20) << "SYSTEM CPU TIME: " << system_cpu() << " sec\n\n";

  if (num_intervals > 1) {
    out << "OVER " << num_intervals << " INTERVALS: " << '\n'
        << setw(20) << "WALL CLOCK TIME: " << sum_wall_clock() << " sec.\n"
        << setw(20) << "  USER CPU TIME: " << sum_user_cpu() << " sec." << '\n'
        << setw(20) << "SYSTEM CPU TIME: " << sum_system_cpu() << " sec.\n\n";
  }
  out.flush();
}

//------------------------------------------------------------------------------------------------//
//! Print out a timing report as a single line summary.
void Timer::printline(std::ostream &out, unsigned const p, unsigned const w) const {
  using std::ios;
  using std::setw;

  out.setf(ios::fixed, ios::floatfield);
  out.precision(p);

  // Width of first column (intervals) should be set by client before calling this function.
  out << num_intervals << setw(w) << sum_user_cpu() << setw(w) << sum_system_cpu() << setw(w)
      << sum_wall_clock() << std::endl;
  out.flush();
}

//------------------------------------------------------------------------------------------------//
// Is this an MPI or Posix timer?
//------------------------------------------------------------------------------------------------//
bool Timer::setIsMPIWtimeAvailable() const {
#ifdef C4_SCALAR
  return false;
#else
  return true;
#endif
}

//------------------------------------------------------------------------------------------------//
// Statics
//------------------------------------------------------------------------------------------------//

/* static */
void Timer::initialize(int & /*argc*/, char ** /*argv*/) {}

//------------------------------------------------------------------------------------------------//
//! Wait until the wall_clock value exceeds the requested pause time.
void Timer::pause(double const pauseSeconds) {
  Require(pauseSeconds > 0.0);

  //! POSIX tms structure for beginning time.
  DRACO_TIME_TYPE tms_begin{};
  //! POSIX tms structure for ending time.
  DRACO_TIME_TYPE tms_end{};

  double begin = wall_clock_time(tms_begin);
  double elapsed(0);

  while (elapsed < pauseSeconds) {
    elapsed = wall_clock_time(tms_end) - begin;
  }
  Ensure(elapsed >= pauseSeconds);
  return;
}

//------------------------------------------------------------------------------------------------//
/*! Print out a summary timing report for averages across MPI ranks.
 *
 * \param out Stream to which to write the report.
 * \param p Precision with which to write the timing and variance numbers.  Defaults to 2.
 * \param w Width of the timing number fields. Defaults to each field being 13 characters wide.
 * \param v Width of the variance number fields. Defaults to each field being 5 characters wide.
 */
void Timer::printline_mean(std::ostream &out, unsigned const p, unsigned const w,
                           unsigned const v) const {
  using std::ios;
  using std::setw;

  unsigned const ranks = rtt_c4::nodes();
  double ni = num_intervals;
  double ni2 = ni * ni;
  double u = sum_user_cpu();
  double u2 = u * u;
  double s = sum_system_cpu();
  double s2 = s * s;
  double ww = sum_wall_clock();
  double ww2 = ww * ww;

  std::array<double, 8> buffer = {ni, ni2, u, u2, s, s2, ww, ww2};
  rtt_c4::global_sum(buffer.data(), 8);

  ni = buffer[0];
  ni2 = buffer[1];
  u = buffer[2];
  u2 = buffer[3];
  s = buffer[4];
  s2 = buffer[5];
  ww = buffer[6];
  ww2 = buffer[7];

  // Casting from a double to unsigned. Ensure that we aren't overflowing the unsigned or dropping a
  // negative sign.
  Check(ni >= 0.0);
  Check(ni < ranks * std::numeric_limits<unsigned>::max());
  auto mni = static_cast<unsigned>(ni / ranks);
  double mu = u / ranks;
  double ms = s / ranks;
  double mww = ww / ranks;

  if (rtt_c4::node() == 0) {
    out.setf(ios::fixed, ios::floatfield);
    out.precision(p);

    // Width of first column (intervals) should be set by client before calling this function.
    out << setw(w) << mni << " +/- " << setw(v)
        << sqrt(fabs(ni2 - 2 * mni * ni + ranks * mni * mni) / ranks) << setw(w) << mu << " +/- "
        << setw(v) << sqrt(fabs(u2 - 2 * mu * u + ranks * mu * mu) / ranks) << setw(w) << ms
        << " +/- " << setw(v) << sqrt(fabs(s2 - 2 * ms * s + ranks * ms * ms) / ranks) << setw(w)
        << mww << " +/- " << setw(v) << sqrt(fabs(ww2 - 2 * mww * ww + ranks * mww * mww) / ranks);
    out << std::endl;
  }
}

} // end namespace rtt_c4

//------------------------------------------------------------------------------------------------//
// end of Timer.cc
//------------------------------------------------------------------------------------------------//
