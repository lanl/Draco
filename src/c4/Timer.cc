//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   c4/Timer.cc
 * \author Thomas M. Evans
 * \date   Mon Mar 25 17:56:11 2002
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "Timer.hh"
#include "ds++/XGetopt.hh"
#include <array>
#include <cmath>
#include <cstdlib>
#include <iomanip>

namespace rtt_c4 {

#ifdef HAVE_PAPI

//------------------------------------------------------------------------------------------------//
/*!
 * Initialize global static non-const data members found in the Timer class
 *
 * By default, we count up the total L2 data cache misses and hits and the total number of floating
 * point operations. These allow us to report the percentage of data cache hits and the number of
 * floating point operations per cache miss.
 */
int Timer::papi_events_[papi_max_counters_] = {PAPI_L2_DCM, PAPI_L2_DCH, PAPI_FP_OPS};
unsigned Timer::papi_num_counters_ = 0;
long long Timer::papi_raw_counts_[papi_max_counters_] = {0, 0, 0};
int selected_cache = 2;

#endif

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

#ifdef HAVE_PAPI

  // Initialize the PAPI library on construction of first timer if it has not already be initialized
  // through a call to Timer::initialize.

  papi_init_();

  // sum of papi wall clock cycles
  sum_papi_wc_cycle = 0;
  // sum of papi wall clock time (microseconds)
  sum_papi_wc_usec = 0;
  // sum of papi virtual cycles
  sum_papi_virt_cycle = 0;
  // sum of papi virtual time (microseconds)
  sum_papi_virt_usec = 0;
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

#ifdef HAVE_PAPI
  double const miss = sum_cache_misses();
  double const hit = sum_cache_hits();
  out << "PAPI Events:\n"
      << setw(26) << 'L' << selected_cache << " cache misses  : " << sum_cache_misses() << "\n"
      << setw(26) << 'L' << selected_cache << " cache hits    : " << sum_cache_hits() << "\n"
      << setw(26) << "Percent hit      : " << 100.0 * hit / (miss + hit) << "\n"
      << setw(26) << "FP operations    : " << sum_floating_operations() << "\n"
      << setw(26) << "Wall Clock cycles: " << sum_papi_wc_cycles() << "\n"
      << setw(26) << "Wall Clock time (us): " << sum_papi_wc_usecs() << "\n"
      << setw(26) << "Virtual cycles: " << sum_papi_virt_cycles() << "\n"
      << setw(26) << "Virtual time (us): " << sum_papi_virt_usecs() << "\n"
      << std::endl;
#endif

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
      << sum_wall_clock();

#ifdef HAVE_PAPI
  double const miss = sum_cache_misses();
  double const hit = sum_cache_hits();
  out << setw(w) << 100.0 * hit / (miss + hit);
  if (papi_num_counters_ > 2) {
    out << setw(w) << sum_floating_operations() / miss;
  }
#endif

  out << std::endl;
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
#ifdef HAVE_PAPI
void Timer::initialize(int &argc, char *argv[])
#else
void Timer::initialize(int & /*argc*/, char ** /*argv*/)
#endif
{
// The initialize function need not be called if the default settings are okay. Otherwise,
// initialize is called with the command line arguments to allow command line control of which cache
// is sampled under PAPI.
//
// At present, there are no non-PAPI options controlled by initialize.
#ifdef HAVE_PAPI
  int j = 0;

  // rtt_dsxx::optind=1; // resets global counter (see XGetopt.cc)

  std::map<std::string, char> long_option;
  long_option["cache"] = 'c';

  int c(0);
  while ((c = rtt_dsxx::xgetopt(argc, argv, (char *)"c:", long_option)) != -1) {
    switch (c) {
    case 'c': // --cache
      char *endptr;
      selected_cache = strtol(argv[i + 1], &endptr, 10);
      if (*endptr != '\0' || selected_cache < 1 || selected_cache > 3) {
        throw std::invalid_argument(" --cache selection is not 1, 2, or 3");
      } else {
        i++;
        switch (selected_cache) {
        case 1:
          papi_events_[0] = PAPI_L1_DCM;
          papi_events_[1] = PAPI_L1_DCH;
          break;

        case 2:
          /* default; no action needed */
          break;

        case 3:
          papi_events_[0] = PAPI_L3_DCM;
          papi_events_[1] = PAPI_L3_DCH;
          break;
        }
      }
      else {
        if (j != i) {
          argv[j] = argv[i];
        }
        j++;
        break;
      }
    default:
      break; // nothing to do.
    }
  }

  int orig_argc = argc;
  argc = j;
  for (; j < orig_argc; ++j) {
    argv[j] = NULL;
  }
#endif // HAVE_PAPI
}

#ifdef HAVE_PAPI
//------------------------------------------------------------------------------------------------//
/* static */
void Timer::papi_init_() {
  static bool first_time = true;
  if (first_time) {
    int retval = PAPI_library_init(PAPI_VER_CURRENT);
    if (retval != PAPI_VER_CURRENT) {
      std::cout << "PAPI library init error!" << std::endl;
      exit(EXIT_FAILURE);
    }
    first_time = false;
  } else {
    int result = PAPI_stop_counters(papi_raw_counts_, papi_num_counters_);
    if (result != PAPI_OK) {
      std::cout << "Failed to stop hardware counters with error " << result << std::endl;

      exit(EXIT_FAILURE);
    }
  }

  if (PAPI_query_event(PAPI_FP_OPS) != PAPI_OK) {
    std::cout << "PAPI: No floating operations counter" << std::endl;
  }

  if (PAPI_query_event(papi_events_[0]) != PAPI_OK) {
    std::cout << "PAPI: No cache miss counter" << std::endl;
  }

  if (PAPI_query_event(papi_events_[1]) != PAPI_OK) {
    std::cout << "PAPI: No cache hit counter" << std::endl;
  }

  papi_num_counters_ = PAPI_num_counters();
  if (papi_num_counters_ < sizeof(papi_events_) / sizeof(int)) {
    std::cout << "PAPI: This system has only " << papi_num_counters_ << " hardware counters.\n"
              << std::endl;
    std::cout << "Some performance statistics will not be available." << std::endl;
  }

  // At present, some platforms *lie* about how many counters they have available, reporting they
  // have three then returning an out of counters error when you actually try to assign the three
  // counter types listed above. Until we have a fix, hardwire to leave out the flops count, which
  // is the least essential of the three counts.
  papi_num_counters_ = 2;

  if (papi_num_counters_ > sizeof(papi_events_) / sizeof(int))
    papi_num_counters_ = sizeof(papi_events_) / sizeof(int);

  int result = PAPI_start_counters(papi_events_, papi_num_counters_);
  if (result != PAPI_OK) {
    std::cout << "Failed to start hardware counters with error " << result << std::endl;

    exit(EXIT_FAILURE);
  }
}
#endif // HAVE_PAPI

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

    // Omit PAPI for now.

    out << std::endl;
  }
}

} // end namespace rtt_c4

//------------------------------------------------------------------------------------------------//
// end of Timer.cc
//------------------------------------------------------------------------------------------------//
