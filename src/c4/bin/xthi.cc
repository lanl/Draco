//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   c4/bin/xthi.cc
 * \author Mike Berry <mrberry@lanl.gov>, Kelly Thompson <kgt@lanl.gov>
 * \date   Wednesday, Aug 09, 2017, 11:45 am
 * \brief  Print MPI rank, thread number and core affinity bindings.
 * \note   Copyright (C) 2017-2025 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "dsxx/SystemCall.hh"
#include "c4/C4_Functions.hh"
#include "c4/xthi_cpuset.hh"
#include <iomanip>
#include <iostream>
#include <string>

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {

  rtt_c4::initialize(argc, argv);
  int const rank = rtt_c4::node();
  std::string const hostname = rtt_dsxx::draco_gethostname();
  unsigned const num_cpus = omp_get_num_procs(); // NOLINT(misc-include-cleaner)

// KPL: See https://gcc.gnu.org/gcc-9/porting_to.html#ompdatasharing for info on the following
// version-specific gcc macro.
#if defined(__GNUC__) && !defined(__clang__) && __GNUC__ < 9
#pragma omp parallel default(none) shared(std::cout)

#else // gcc version >= 9
#pragma omp parallel default(none) shared(num_cpus, hostname, rank, std::cout)
#endif
  {
    int const thread = omp_get_thread_num(); // NOLINT(misc-include-cleaner)
    std::string const cpuset = rtt_c4::cpuset_to_string(num_cpus);

#pragma omp critical
    {
      std::cout << hostname << " :: Rank " << std::setfill('0') << std::setw(5) << rank
                << ", Thread " << std::setfill('0') << std::setw(3) << thread
                << ", core affinity = " << cpuset << '\n';
    } // end omp critical
  }   // end omp parallel

  rtt_c4::finalize();
  return (0);
}

//------------------------------------------------------------------------------------------------//
// End c4/bin/xthi.cc
//------------------------------------------------------------------------------------------------//
