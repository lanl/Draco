//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   diagnostics/Timing.hh
 * \author T.M. Kelly, Thomas M. Evans
 * \date   Tue Dec 13 10:44:29 2005
 * \brief  Timing class and macros definition.
 * \note   Copyright (C) 2010-2021 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef diagnostics_Timing_hh
#define diagnostics_Timing_hh

#include "diagnostics/config.h"
#include "ds++/config.h"
#include <map>
#include <string>
#include <vector>

namespace rtt_diagnostics {

//================================================================================================//
/*!
 * \class Timing_Diagnostics
 * \brief Class to hold timing results for diagnostic output.
 *
 * This class provides a simple interface to store timing data during a simulation.  Generally, one
 * adds a timer label and value using the update_timer() function:
 * \code
 *   // start a timer
 *   ...
 *   // do work
 *   ...
 *   // stop timer
 *   rtt_diagnostics::Timing_Diagnostics::update_timer("Solver", time);
 * \endcode
 *
 * There is no need to "add" the timer entry for "Solver" before an update.  If the key "Solver"
 * does not exist it is added with value 0.0 before applying the value.
 *
 * The easiest way to use this class is through the TIMER macros.
 *
 * \sa \ref diagnostics_timing
 *
 * \example diagnostics/test/tstTiming.cc
 */
//================================================================================================//

class Timing_Diagnostics {
public:
  // Useful typedef.
  using Vec_Keys = std::vector<std::string>;

private:
  // >>> PRIVATE DATA MEMBERS

  //! Map of timers.
  DLL_PUBLIC_diagnostics static std::map<std::string, double> timers;

public:
  // >>> FUNCTIONAL INTERFACE

  // Add a value to the timer with name key.
  static void update_timer(const std::string &key, double value);

  //! Get a timer's value.  Adds timer with name key to map.
  static double timer_value(const std::string &k) { return timers[k]; }

  //! Get number of timers in map.
  static size_t num_timers() { return timers.size(); }

  // Return a vector of timer keys.
  static Vec_Keys timer_keys();

  // Reset a timer.
  static void reset_timer(const std::string &key);

  // Reset all timers from the map of timers.
  static void reset_timers();

  // Delete a timer from the map of timers.
  static void delete_timer(const std::string &key);

  // Delete all timers from the map of timers.
  static void delete_timers();

  // >>> IMPLEMENTATION

  // Disable default ctor
  Timing_Diagnostics() = delete;

  // Disable default dtor
  ~Timing_Diagnostics() = delete;
};

} // end namespace rtt_diagnostics

//------------------------------------------------------------------------------------------------//
/*!
 * \page diagnostics_timing Macros for timing
 *
 * Four macros are defined here; these macros insert timer calls into code if a global definition,
 * \c DRACO_TIMING, is greater then 0.  The default value is to set DRACO_TIMING == 0. They use the
 * draco rtt_c4::Timer and rtt_diagnostics::Timing_Diagnostics classes.
 *
 * As of June 2020, the macros can be implemented in two different ways. One may use the Caliper
 * library if that is found through the Draco build system; or one may continue to use the Draco
 * timing diagnostics. The Caliper approach will probably be more useful in the long term.
 *
 * The build system sets DRACO_TIMING. The following settings apply:
 * - 0 turns off all TIMER macros
 * - 1 turns on TIMER, TIMER_START, TIMER_STOP, and TIMER_RECORD
 * - 2 turns on all TIMER macros (include TIMER_REPORT)
 * .
 * The default is 0.
 *
 * In code,
 * \code
 * #include "diagnostics/Timing.hh"
 *
 * TIMER(foo);
 * TIMER_START("Snippet", foo);
 * // ...
 * // code interval to time
 * // ...
 * TIMER_STOP("Snippet", foo);
 * TIMER_RECORD( "Snippet", foo);
 * TIMER_REPORT( foo, std::cout, "interval 42");
 * \endcode
 * would produce the following output:
 * \verbatim
 *   somefile.cc ###: interval 42 elapsed wall_clock: ## seconds; elapsed user_time: ## seconds;
 *   elapsed sys_time ## seconds.
 * \endverbatim
 * The key "Snippet" can be used to access the stored time through the Timer_Diagnostics class:
 * \code
 * #ifdef DRACO_TIMING_ON
 *   vector<string> keys = Timing_Diagnostics::timer_keys();
 *   for (int i = 0; i < keys.size(); i++)
 *       cout << keys[i] << "\t" << Timing_Diagnostics::timer_value(keys[i]) << endl;
 *   Timing_Diagnostics::reset_timers();
 * #endif
 * \endcode
 *
 * \dot
 * digraph draco_timing {
 *  node [shape=record, fontname=Helvetica, fontsize=10];
 *      start [ label="Start", shape="box", style=rounded ];
 *      dt_gt_0 [ label="DRACO_TIMING>0?", shape="diamond", URL="\ref DRACO_TIMING"];
 *      have_caliper [ label="DRACO_CALIPER?", shape="diamond", URL="\ref DRACO_CALIPER" ];
 *      notimers [ label="No timers", shape="box", style=rounded ];
 *      caltimer [ label="Caliper Timers", shape="box", style=rounded ];
 *      dratimer [ label="Draco Timers", shape="box", style=rounded ];
 *
 *      start -> dt_gt_0 [ arrowhead="open", style="dashed" ];
 *      dt_gt_0 -> notimers [label="no"];
 *      dt_gt_0 -> have_caliper [label="yes"];
 *      have_caliper -> caltimer [label="yes"];
 *      have_caliper -> dratimer [label="no"];
 * }
 * \enddot
 */

//------------------------------------------------------------------------------------------------//
/*!
 * \def TIMER(timer_name)
 *
 * If DRACO_TIMING_ON is defined, TIMER( timer_name) expands to:
 * \code
 *     rtt_c4::Timer timer_name
 * \endcode
 * Otherwise it is empty.
 */

//------------------------------------------------------------------------------------------------//
/*!
 * \def TIMER_START(segment_name, timer_name)
 *
 * If DRACO_TIMING > 0 and DRACO_CALIPER is false. TIMER_START(segment_name, timer_name) expands to:
 * \code
 *     timer_name.start()
 * \endcode
 *
 * (Note that the segment_name is ignored.) If DRACO_TIMING > 0 and DRACO_CALIPER is true, then
 * TIMER_START(segment_name, timer_name) expands to:
 *
 * \code
 *     CALI_MARK_BEGIN(segment_name)
 * \endcode
 * (Note that the timer_name is ignored.). Otherwise the macro expansion is empty.
 */

//------------------------------------------------------------------------------------------------//
/*!
 * \def TIMER_STOP(segment_name, timer_name)
 *
 * If DRACO_TIMING_ON > 0, and DRACO_CALIPER is false, then TIMER_STOP(segment_name, timer_name)
 * expands to:
 *
 * \code
 *     timer_name.stop()
 * \endcode
 *
 * (Note that the segment_name is ignored.) If DRACO_TIMING > 0 and DRACO_CALIPER is true, then
 * TIMER_STOP(segment_name, timer_name) expands to:
 *
 * \code
 *     CALI_MARK_END(segment_name)
 * \endcode
 *
 * (Note that the timer_name is ignored.). Otherwise the macro expansion is empty. Otherwise it is
 * empty.
 */

//------------------------------------------------------------------------------------------------//
/*!
 * \def TIMER_RECORD( name, timer)
 *
 * If DRACO_TIMING_ON > 0, and DRACO_CALIPER is false TIMER_RECORD( name, timer) expands to:
 *
 * \code
 *     rtt_diagnostics::Timing_Diagnostics::update_timer(name, timer.wall_clock())
 * \endcode
 *
 * Otherwise it is empty.
 */

//------------------------------------------------------------------------------------------------//
/*!
 * \def TIMER_REPORT( timer_name, ostream, comment)
 *
 * If DRACO_TIMING > 1, , and DRACO_CALIPER is false TIMER_REPORT( timer_name, ostream, comment)
 * expands to:
 *
 * \code
 *     ostream << __FILE__ << " " << __LINE__ << ": " << comment      \
 *             << " elapsed wall_clock: " << timer.wall_clock()       \
 *             << " seconds; elapsed user_time: " << timer.user_cpu() \
 *             << " seconds; elapsed sys_time: " << timer.system_cpu()\
 *             << " seconds.\n" << flush
 * \endcode
 *
 * Otherwise it is empty. The flush ensures that regression tests continue to pass (otherwise, in
 * parallel runs, output may arrive "out of order" and trample the output that the regression tests
 * look for).
 *
 * \sa \ref diagnostics_timing
 *
 * Activation
 *
 * \code
 * if DRACO_TIMING == 0
 *   // all timers inactive
 * else
 *   // all timers are available
 *   ifdef DRACO_CALIPER
 *     // use Caliper timers
 *   else
 *     // use Draco native timers
 *   endif
 * endif
 * \endcode
 */
//------------------------------------------------------------------------------------------------//

#if !defined(DRACO_TIMING)
#define DRACO_TIMING 0
#endif

//------------------------------------------------------------------------------------------------//
#if DRACO_TIMING == 0

//! All timing operations are inactive.

#define TIMER(timer)
#define TIMER_START(name, timer)
#define TIMER_STOP(name, timer)
#define TIMER_RECORD(name, timer)
#define TIMER_REPORT(timer, ostream, comment)

#endif

//------------------------------------------------------------------------------------------------//
//! DRACO NATIVE TIMERS (not caliper)
//------------------------------------------------------------------------------------------------//

#ifndef DRACO_CALIPER
#if DRACO_TIMING > 0

#include "c4/Timer.hh"

#define DRACO_TIMING_ON
#define TIMER(timer) rtt_c4::Timer timer
#define TIMER_START(name, timer) timer.start()
#define TIMER_STOP(name, timer) timer.stop()
#define TIMER_RECORD(name, timer)                                                                  \
  rtt_diagnostics::Timing_Diagnostics::update_timer(name, timer.wall_clock())

#endif // DRACO_TIMING > 0

#if DRACO_TIMING > 1

//! Turn on timing report output.  This is an add-on option to the basic timing operations.
#define TIMER_REPORT(timer, ostream, comment)                                                      \
  ostream << __FILE__ << " " << __LINE__ << ": " << comment                                        \
          << " elapsed wall_clock: " << timer.wall_clock()                                         \
          << " seconds; elapsed user_time: " << timer.user_cpu()                                   \
          << " seconds; elapsed sys_time: " << timer.system_cpu() << " seconds.\n"                 \
          << std::flush

#else // DRACO_TIMING > 1
#define TIMER_REPORT(timer, ostream, comment)
#endif // DRACO_TIMING > 1

//------------------------------------------------------------------------------------------------//
// CALIPER TIMERS
//------------------------------------------------------------------------------------------------//

#else // Caliper is available
#include <caliper/cali.h>

#if DRACO_TIMING > 0

#include "c4/Timer.hh"
#define DRACO_TIMING_ON
#define TIMER(timer)
#define TIMER_START(name, timer) CALI_MARK_BEGIN(name)
#define TIMER_STOP(name, timer) CALI_MARK_END(name)
#define TIMER_RECORD(name, timer)
#define TIMER_REPORT(timer, ostream, comment)

#endif // DRACO_TIMING > 0

#endif // DRACO_CALIPER

#endif // diagnostics_Timing_hh

//------------------------------------------------------------------------------------------------//
//  end of diagnostics/Timing.hh
//------------------------------------------------------------------------------------------------//
