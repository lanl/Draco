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
 *      start [ label="Start", shape="box", style="rounded" ];
 *      dt_gt_0 [ label="DRACO_TIMING>0?", shape="diamond", URL="\ref DRACO_TIMING"];
 *      have_caliper [ label="DRACO_CALIPER?", shape="diamond" ];
 *      notimers [ label="No timers", shape="box", style="rounded" ];
 *      caltimer [ label="Caliper Timers", shape="box", style="rounded" ];
 *      dratimer [ label="Draco Timers", shape="box", style="rounded" ];
 *
 *      start -> dt_gt_0 [ arrowhead="open", style="dashed" ];
 *      dt_gt_0 -> notimers [label="no"];
 *      dt_gt_0 -> have_caliper [label="yes"];
 *      have_caliper -> caltimer [label="yes"];
 *      have_caliper -> dratimer [label="no"];
 * }
 * \enddot
 *
 * If Caliper timers are used, simple reports can be generated like this:
 *
 * \verbatim
 $ CALI_CONFIG_PROFILE=runtime-report ./tstTiming

  Path            Inclusive time Exclusive time Time %
  Outer                 0.000710       0.000014  1.680672
    C_iteration         0.000116       0.000116 13.925570
    B_iteration         0.000231       0.000231 27.731092
    A_iteration         0.000349       0.000349 41.896759
 * \endverbatim
 *
 * Or with MPI like this:
 *
 * \verbatim
$ CALI_CONFIG_PROFILE=runtime-report mpiexec -n 2 ../bin/mcgrid -i amr_xyz_mg/mginf01
...
Path                                 Inclusive time Exclusive time Time %
execute_IMC:regroup                        0.006950       0.006950  9.675084
execute_IMC:run_IMC                        0.012787       0.003883  5.405518
  host_loop                                0.008904       0.000609  0.847788
    trasport_time                          0.008295       0.008295 11.547457
execute_IMC:cycle_init                     0.024644       0.005712  7.951666
  Build_Source                             0.011629       0.008668 12.066709
    build_Source:comb_census               0.000798       0.000027  0.037587
      rep_comb_census                      0.000771       0.000771  1.073308
    build_Source:calc_src_num              0.000443       0.000443  0.616700
    build_Source:build_strategy            0.000593       0.000474  0.659855
      calc_evol                            0.000119       0.000119  0.165660
    build_Source:init_cen_nrg              0.001127       0.001087  1.513211
      calc_evol                            0.000040       0.000040  0.055684
  Build_Mat_State                          0.004806       0.004806  6.690425
  cycle_init:fine_tune_census_cell         0.002497       0.002497  3.476070
 * \endverbatim
 *
 * You can also generate MPI timing reports like this:
 *
 * \verbatim
 $ CALI_CONFIG_PROFILE=mpi-runtime-report mpiexec -n 2 ../bin/mcgrid -i amr_xyz_mg/mginf01
...
Path                                 Min time/rank Max time/rank Avg time/rank Time % (total)
MPI_Finalize                              0.000007      0.000042      0.000024       0.036926
MPI_Comm_free                             0.000048      0.000054      0.000051       0.076866
MPI_Type_size                             0.000006      0.000007      0.000006       0.009797
MPI_Type_commit                           0.000008      0.000009      0.000008       0.012811
MPI_Type_create_struct                    0.000022      0.000023      0.000022       0.033911
MPI_Scan                                  0.000013      0.000022      0.000018       0.026375
MPI_Allreduce                             0.000083      0.016863      0.008473      12.770254
MPI_Barrier                               0.000091      0.006369      0.003230       4.868160
MPI_Type_free                             0.000017      0.000018      0.000017       0.026375
execute_IMC:regroup                       0.001608      0.001620      0.001614       2.432573
  MPI_Allreduce                           0.000395      0.000408      0.000402       0.605129
execute_IMC:run_IMC                       0.003570      0.003671      0.003621       5.456710
  host_loop                               0.000556      0.000559      0.000557       0.840247
    trasport_time                         0.007417      0.007493      0.007455      11.235955
execute_IMC:cycle_init                    0.005869      0.005936      0.005903       8.896073
  Build_Source                            0.002678      0.002696      0.002687       4.049767
    MPI_Send                              0.000029      0.000029      0.000029       0.021854
    build_Source:comb_census              0.000020      0.000020      0.000020       0.030143
      MPI_Allreduce                       0.000010      0.000012      0.000011       0.016579
      rep_comb_census                     0.000618      0.000657      0.000638       0.960821
        MPI_Allreduce                     0.000028      0.000066      0.000047       0.070837
    build_Source:calc_src_num             0.000384      0.000402      0.000393       0.592318
      MPI_Send                            0.000020      0.000020      0.000020       0.015072
      MPI_Barrier                         0.000023      0.000065      0.000044       0.066315
      MPI_Get_count                       0.000013      0.000013      0.000013       0.009797
      MPI_Recv                            0.000030      0.000030      0.000030       0.022608
    build_Source:build_strategy           0.000406      0.000420      0.000413       0.622461
      MPI_Send                            0.000012      0.000012      0.000012       0.009043
      MPI_Barrier                         0.000015      0.000043      0.000029       0.043708
      MPI_Get_count                       0.000010      0.000010      0.000010       0.007536
      MPI_Recv                            0.000025      0.000025      0.000025       0.018840
      MPI_Allreduce                       0.000039      0.000050      0.000044       0.067069
      calc_evol                           0.000093      0.000096      0.000094       0.142428
    MPI_Allreduce                         0.000302      0.000306      0.000304       0.458180
    build_Source:init_cen_nrg             0.000544      0.000553      0.000549       0.826683
      MPI_Send                            0.000010      0.000010      0.000010       0.007536
      MPI_Barrier                         0.000012      0.000033      0.000022       0.033911
      MPI_Get_count                       0.000007      0.000007      0.000007       0.005275
      MPI_Recv                            0.000021      0.000021      0.000021       0.015825
      calc_evol                           0.000040      0.000042      0.000041       0.061794
    MPI_Barrier                           0.000022      0.000093      0.000058       0.086662
    MPI_Get_count                         0.000010      0.000010      0.000010       0.007536
    MPI_Recv                              0.005526      0.005526      0.005526       4.164312
  Build_Mat_State                         0.004419      0.004634      0.004526       6.822207
  cycle_init:fine_tune_census_cell        0.001510      0.001541      0.001525       2.299188
    MPI_Allreduce                         0.000039      0.000100      0.000069       0.104748
 * \endverbatim
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
 * If DRACO_TIMING > 1, and DRACO_CALIPER is false TIMER_REPORT( timer_name, ostream, comment)
 * expands to:
 *
 * \code
 *   ostream << __FILE__ << " " << __LINE__ << ": " << comment << " elapsed wall_clock: " \
 *           << timer.wall_clock() << " seconds; elapsed user_time: " << timer.user_cpu() \
 *           << " seconds; elapsed sys_time: " << timer.system_cpu() << " seconds.\n" << flush
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
/*!
 * DRACO_CALIPER
 *
 * If defined, this CPP macro indicates that Caliper is available in the current build of Draco. The
 * logic for setting this variable is in diagnostics/CMakeLists.txt and the value is safed to
 * diagnostics/config.h.
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
// DRACO NATIVE TIMERS (not caliper)
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
