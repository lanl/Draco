//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file  ds++/SystemCall.hh
 * \brief Wrapper for system calls. Hide differences between Unix/Windows system calls.
 * \note  Copyright (C) 2012-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_dsxx_SystemCall_hh
#define rtt_dsxx_SystemCall_hh

#include "ds++/config.h"
#include <array>
#include <chrono>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <thread>
#ifdef _MSC_VER
#define _WINSOCKAPI_
#include <WinSock2.h>
#include <Windows.h>
#include <sys/types.h>
#endif
#include <sys/stat.h> // stat (UNIX) or _stat (WIN32)

namespace rtt_dsxx {

//! Character used as path separator.
char constexpr WinDirSep = '\\';
char constexpr UnixDirSep = '/';
#ifdef _MSC_VER
char const dirSep = WinDirSep;
std::string const exeExtension(".exe");
#else
char const dirSep = UnixDirSep;
std::string const exeExtension("");
#endif

//================================================================================================//
// FREE FUNCTIONS
//================================================================================================//

//! Return the local hostname
std::string draco_gethostname();

//! Return the local process id
int draco_getpid();

//! Return the current working directory
std::string draco_getcwd();

//! Return the stat value for a file
class draco_getstat {
private:
  int stat_return_code{0};
#ifdef _MSC_VER
  struct _stat buf;
  bool filefound{false};
  WIN32_FIND_DATA FileInformation{}; // Additional file information
#else
  struct stat buf {};
#endif

public:
  //! constructor
  explicit draco_getstat(std::string const &fqName);
  //! If the call to stat failed, this function will return false.
  bool valid() { return stat_return_code == 0; }
  bool isreg();
  bool isdir();
  int errorCode() { return stat_return_code; }
  /*!
   * \brief Determine if the file has the requested permission bits set.
   * \note The leading zero for the mask is important.
   */
  bool has_permission_bit(int mask = 0777);
};

//! Use Linux realpath to resolve symlinks
std::string draco_getrealpath(std::string const &path);

//! Create a directory
void draco_mkdir(std::string const &path);

/*!
 * \brief Remove file or directory (not recursive)
 *
 * For recursive directory delete, see path.hh's walk_directory_tree and the functor wdtOpRemove.
 */
void draco_remove(std::string const &path);

} // namespace rtt_dsxx

/*!
 * \brief set and unset environment variables.
 */
#ifdef MSVC
#define draco_unsetenv(k) _putenv_s(k, "")
#define draco_setenv(k, v) _putenv_s(k, v)
#else
#define draco_unsetenv(k) unsetenv(k)
#define draco_setenv(k, v) setenv(k, v, 1)
#endif

#ifndef MSVC
//================================================================================================//
/*!
 * \class CPUinfo
 * \brief Read information from /proc/stat and provide information about CPU states.
 *
 * \example c4/test/tstOMP.cc
 * The main() reports the average load and cpu utilization (for Linux).
 */
//================================================================================================//
class CPUinfo {

private:
  //! The cpudata has 10 values (ref: /proc/stat)
  static constexpr size_t cpudata_num_values{10};
  /*!
   * \struct CPUData
   *
   * Data struct for cpu information that will be captured.  The format for /proc/stat is something
   * like this example:
   * <code>
   *    cpu  764699971 57752895 158644540 75321983666 67359813 0 736586 0 0 0
   * </code>
   * A string followed by 10 size_t values.
   *
   * The 10 values represent cpu states:
   * - user (0)
   * - nice (1)
   * - system (2)
   * - idle (3)
   * - iowait (4)
   * - irg (5)
   * - softirq (6)
   * - steal (7)
   * - guest (8)
   * - guest_nice (9)
   * .
   */
  struct CPUData {
    std::string cpu{};
    size_t num_cpus{0};
    std::array<size_t, cpudata_num_values> times{};
  };

  /* Maniulators */

  //! Read /proc/stat to determine the system load.
  CPUData ReadStatsCPU() {
    CPUData cpu_all_entry;
    const std::string label_cpu_all("cpu");
    std::ifstream fileStat("/proc/stat");
    std::string line("");
    size_t num_cpu(0);
    while (std::getline(fileStat, line)) { // Read a line
      // cpu stats line found
      if (!line.compare(0, label_cpu_all.length(), label_cpu_all)) {
        std::istringstream ss(line);
        // read cpu label
        ss >> cpu_all_entry.cpu;
        // read times
        if (cpu_all_entry.cpu == "cpu") {
          for (auto &t : cpu_all_entry.times) {
            ss >> t;
          }
        } else {
          num_cpu++;
        }
      }
      cpu_all_entry.num_cpus = num_cpu;
    }
    return cpu_all_entry;
  }

  /* Private accessors */

  size_t machIdleTime(CPUData const &c) const { return c.times[3] + c.times[4]; }
  size_t machActiveTime(CPUData const &c) const {
    return c.times[0] + c.times[1] + c.times[2] + c.times[5] + c.times[6] + c.times[7] +
           c.times[8] + c.times[9];
  }

  //! Read /proc/loadavg data.
  std::array<float, 3> ReadLoadAvg() {
    std::array<float, 3> retvalue{0.0, 0.0, 0.0};
    std::ifstream fileStat("/proc/loadavg");
    std::string line("");
    while (std::getline(fileStat, line)) { // Read a line
      std::istringstream ss(line);
      for (auto &t : retvalue) {
        ss >> t;
      }
    }
    return retvalue;
  }

public:
  //! Report cpu usage
  void report_cpu_utilization() {
    // initial state
    CPUData entry1 = ReadStatsCPU();
    std::this_thread::sleep_for(std::chrono::milliseconds(500));
    // state after 0.1 sec.
    CPUData entry2 = ReadStatsCPU();

    auto const activeTime = static_cast<float>(machActiveTime(entry2) - machActiveTime(entry1));
    auto const idleTime = static_cast<float>(machIdleTime(entry2) - machIdleTime(entry1));
    float const totalTime = activeTime + idleTime;
    std::cout.width(6);
    std::cout.precision(2);
    std::cout << "Average CPU Utilization      : " << 100.0F * activeTime / totalTime << "%\n"
              << "Number of cores (hw threads) :" << entry1.num_cpus << std::endl;
  }

  //! Report load average
  void report_load_average() {
    // 1, 5, and 15 minute load averages reported by /proc/loadavg.
    std::array<float, 3> loadAvg = ReadLoadAvg();
    std::cout << "Load Average: " << loadAvg[0] << ", " << loadAvg[1] << ", " << loadAvg[2]
              << " (1 min, 5 min, 15 min)" << std::endl;
  }
};
#endif

#endif // rtt_dsxx_SystemCall_hh

//------------------------------------------------------------------------------------------------//
// end of SystemCall.hh
//------------------------------------------------------------------------------------------------//
