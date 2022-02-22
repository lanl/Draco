//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   c4/bin/ythi.hh
 * \author Kelly Thompson <kgt@lanl.gov>, Tim Kelley <tkelley@lanl.gov.
 * \date   Tuesday, Jun 05, 2018, 17:12 pm
 * \brief  Print MPI rank, thread number and core affinity bindings.
 * \note   Copyright (C) 2019-2022 Triad National Security, LLC., All rights reserved.
 *
 * These functions are used by c4/bin/ythi.cc and c4/test/tstQuoWrapper.cc
 */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_c4_bin_ythi_hh
#define rtt_c4_bin_ythi_hh

#include "c4/C4_Functions.hh"
#include "c4/SLURM_Task_Info.hh"
#include "c4/xthi_cpuset.hh"
#include <atomic>
#include <iomanip>
#include <iostream>
#include <thread>
#include <vector>

namespace rtt_c4 {

//------------------------------------------------------------------------------------------------//
/**\brief After atomic bool changes to true, print out some thread info. */
void run_thread(std::atomic<bool> &signal, std::string const &hostname, int const rank,
                size_t const simple_thread_id) {
  while (!signal) {
  }
  unsigned const num_cpu = std::thread::hardware_concurrency();
  std::string cpuset = rtt_c4::cpuset_to_string(num_cpu);
  std::cout << hostname << " :: Rank " << std::setfill('0') << std::setw(5) << rank << ", Thread "
            << std::setfill('0') << std::setw(3) << simple_thread_id
            << ", core affinity = " << cpuset << std::endl;
  return;
}

//------------------------------------------------------------------------------------------------//
void report_bindings(uint32_t const num_workers) {

  std::string const hostname = rtt_dsxx::draco_gethostname();
  std::vector<std::atomic<bool>> signals(num_workers);
  std::vector<std::thread> threads(num_workers);
  uint32_t const num_cpus = std::thread::hardware_concurrency();

  for (size_t i = 0; i < num_workers; ++i) {
    signals[i].store(false);
    threads[i] =
        std::thread(run_thread, std::ref(signals[i]), std::ref(hostname), rtt_c4::rank(), i + 1);
  }
  std::string cpuset = rtt_c4::cpuset_to_string(num_cpus);
  int const host_thread(0);
  std::cout << hostname << " :: Rank " << std::setfill('0') << std::setw(5) << rtt_c4::rank()
            << ", Thread " << std::setfill('0') << std::setw(3) << host_thread
            << ", core affinity = " << cpuset << std::endl;
  for (size_t i = 0; i < num_workers; ++i) {
    signals[i].store(true);
    threads[i].join();
  }

  return;
}
} // namespace rtt_c4

#endif // rtt_c4_bin_ythi_hh

//------------------------------------------------------------------------------------------------//
// End c4/bin/ythi.hh
//------------------------------------------------------------------------------------------------//
