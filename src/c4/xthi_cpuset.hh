//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   c4/xthi_cpuset.hh
 * \author Mike Berry <mrberry@lanl.gov>, Kelly Thompson <kgt@lanl.gov>
 * \date   Wednesday, Aug 09, 2017, 11:45 am
 * \brief  Helper functions to generate string for core affinity.
 * \note   Copyright (C) 2019-2021 Triad National Security, LLC., All rights reserved.
 *
 * These functions are needed by c4's xthi and ythi programs to report human readable thread
 * bindings.  It is also used by the unit test for libquo.
 */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_c4_xthi_cpuset_hh
#define rtt_c4_xthi_cpuset_hh

#include "c4_omp.h"
#include "ds++/Assert.hh"
#include "ds++/SystemCall.hh"
#include <bitset>
#include <sstream>

#ifdef _MSC_VER
#include <processthreadsapi.h> // requries SystemCall.hh to be loaded first.
#endif

/*! \bug Broken on OSX (sched_getaffinity, cpu_set_t, etc. not defined).
 *
 * Possible solution at https://yyshen.github.io/2015/01/18/binding_threads_to_cores_osx.html
 */

namespace rtt_c4 {

//------------------------------------------------------------------------------------------------//
/*!
 * \brief query the OS for an affinity bitmask and translate it to a human
 *        readable form.
 *
 * For Linux-like systems, borrow code from util-linux-2.13-pre7/schedutils/taskset.c.
 *
 * For Windows systems, make use of GetProcessAffinityMask (processthreadsapi.h). Refs:
 *
 * - https://docs.microsoft.com/en-us/windows/desktop/api/winbase/nf-winbase-getprocessaffinitymask
 * - https://stackoverflow.com/questions/10877182/getprocessaffinitymask-returns-processaffinty-and-systemaffinity-as-1-overflow
 * - https://stackoverflow.com/questions/2215063/how-can-you-find-the-processor-number-a-thread-is-running-on
 * .
 *
 * In both cases, the affinity bitmask is restricted to 64 cores.  If a node has more than 64 cores,
 * extra logic will be needed to make this function work.
 *
 * \return A string of the form "0-8;16-32;" or "0-63"
 */
#ifdef APPLE

#include <sys/sysctl.h>
#include <sys/types.h>

/*
 * I think this is the mac equivelent of CPU_SETSIZE
 *
 * example - https://stackoverflow.com/questions/150355/programmatically-find-the-number-of-cores-on-a-machine
 * documentation - https://developer.apple.com/library/archive/documentation/System/Conceptual/ManPages_iPhoneOS/man3/sysconf.3.html
 *
*/
#ifndef CPU_SETSIZE
#define CPU_SETSIZE _SC_NPROCESSORS_ONLN
#endif

/*! \bug Broken on OSX (sched_getaffinity, cpu_set_t, etc. not defined).
 *
 * Possible solution at https://yyshen.github.io/2015/01/18/binding_threads_to_cores_osx.html
 */

#define SYSCTL_CORE_COUNT "machdep.cpu.core_count"

using cpu_set_t = struct cpu_set { uint64_t count; };

static inline void CPU_ZERO(cpu_set_t *cs) { cs->count = 0; }

static inline void CPU_SET(int num, cpu_set_t *cs) { cs->count |= (1 << num); }

static inline bool CPU_ISSET(int num, cpu_set_t *cs) {
  uint64_t teh_bit = 1ull << num;
  return (cs->count & teh_bit) != 0;
}

inline int sched_getaffinity(pid_t /*pid*/, size_t /*cpu_size*/, cpu_set_t *cpu_set) {
  int64_t core_count = 0;
  size_t len = sizeof(core_count);
  int ret = sysctlbyname(SYSCTL_CORE_COUNT, &core_count, &len, 0, 0);
  if (ret) {
    printf("error while get core count %d\n", ret);
    return -1;
  }
  cpu_set->count = 0;
  for (int i = 0; i < core_count; i++) {
    cpu_set->count |= (1 << i);
  }

  return 0;
}

#endif

#ifdef _MSC_VER

//! \param[in] num_cpu Number of CPU's per node.
inline std::string cpuset_to_string(unsigned const num_cpu) {

  // return value;
  std::ostringstream cpuset;
  // The thread affinity bitmask functions used below are limited to 64.
  Insist(num_cpu <= 64, "Might need to use alternate cpu-groups information with this fuction!");

  DWORD_PTR dwProcessAffinity;
  DWORD_PTR dwSystemAffinity;
  bool ok = GetProcessAffinityMask(GetCurrentProcess(), &dwProcessAffinity, &dwSystemAffinity);
  Insist(ok, "GetProcessAffinityMask() failed!");

  // Convert the bitmask to an array of bools and then to a string to represent a CPU range
  std::bitset<64> const affmask(dwProcessAffinity);
  size_t begin(0);
  bool enabled(false);
  for (size_t i = 0; i < affmask.size(); ++i) {
    if (enabled) {
      if (!affmask[i]) {
        enabled = false;
        cpuset << begin << "-" << i - 1 << "; ";
      }
    } else {
      // looking for next core that can be used
      if (affmask[i]) {
        enabled = true;
        begin = i;
      }
    }
  }
  return cpuset.str();
}

#else

inline std::string cpuset_to_string(unsigned const /*num_cpu*/) {

  // return value;
  std::ostringstream cpuset;
  // local storage; retrieve the thread affinity bitmask
  cpu_set_t coremask;
  (void)sched_getaffinity(0, sizeof(coremask), &coremask);

  // Convert the bitmask into something that is human readable.
  size_t entry_made = 0;

  for (int i = 0; i < CPU_SETSIZE; i++) {
    if (CPU_ISSET(i, &coremask)) {
      int run = 0;
      entry_made = 1;
      for (int j = i + 1; j < CPU_SETSIZE; j++) {
        if (CPU_ISSET(j, &coremask))
          run++;
        else
          break;
      }
      if (run == 0)
        cpuset << i << ",";
      else
        cpuset << i << "-" << i + run << ",";
      i += run;
    }
  }
  return cpuset.str().substr(0, cpuset.str().length() - entry_made);
}

#endif

} // end namespace rtt_c4

#endif // rtt_c4_xthi_cpuset_hh

//------------------------------------------------------------------------------------------------//
// End c4/xthi_cpuset.hh
//------------------------------------------------------------------------------------------------//
