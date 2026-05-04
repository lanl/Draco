//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   c4/bin/ythi.cc
 * \author Mike Berry <mrberry@lanl.gov>, Kelly Thompson <kgt@lanl.gov>,
 *         Tim Kelley <tkelley@lanl.gov.
 * \brief  Print MPI rank, thread number and core affinity bindings.
 * \note   Copyright (C) 2018-2026 Triad National Security, LLC., All rights reserved.
 *
 * Rewritten by Tim Kelley to run C++11 std::threads You may override \c NUM_WORKERS on the compile
 * command line.  For example to run with 4 worker threads:
 *
 * \code
 * $ ./ythi 4
 * \endcode
 *
 * The default is 1 worker thread (over and above the host thread)
 *
 * If executed from a SLURM allocation and the \c NUM_WORKERS value is not set, then the number of
 * workers will be set to the environment variable `SLURM_CPUS_PER_TASK`, which is set by SLURM.
 */
//------------------------------------------------------------------------------------------------//

// cSpell:ignore SLURM, ythi

#include "c4/bin/ythi.hh"
#include "c4/C4_Functions.hh"
#include "c4/SLURM_Task_Info.hh"
#include <cstdint>
#include <iostream>
#include <string>
#include <thread>

//------------------------------------------------------------------------------------------------//
int main(int argc, char **argv) {
  rtt_c4::SLURM_Task_Info const sti;
  uint32_t const numthreads = sti.is_cpus_per_task_set() ? sti.get_cpus_per_task() : 1;
  uint32_t const YTHI_NUM_WORKERS = (argc > 1) ? std::stoi(argv[1]) : numthreads - 1;
  unsigned const num_cpus = std::thread::hardware_concurrency();

  rtt_c4::initialize(argc, argv);
  {
    uint32_t const myRank = rtt_c4::rank();
    if (myRank == 0)
      std::cout << "Found " << num_cpus << " logical CPUs (hardware therads) per node.\n";
    rtt_c4::report_bindings(YTHI_NUM_WORKERS);
  }
  rtt_c4::finalize();
  return (0);
}

//------------------------------------------------------------------------------------------------//
// End c4/bin/ythi.cc
//------------------------------------------------------------------------------------------------//
