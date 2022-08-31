//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   c4/bin/xthi.cc
 * \author Mike Berry <mrberry@lanl.gov>, Kelly Thompson <kgt@lanl.gov>
 * \date   Wednesday, Aug 09, 2017, 11:45 am
 * \brief  Print MPI rank, thread number and core affinity bindings.
 * \note   Copyright (C) 2017-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "c4/C4_Functions.hh"
#include "c4/xthi_cpuset.hh"
#include <iomanip>
#include <iostream>

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {

  rtt_c4::initialize(argc, argv);
  int const rank = rtt_c4::node();
  std::string const hostname = rtt_dsxx::draco_gethostname();
  unsigned const num_cpus = omp_get_num_procs();

// KPL: See https://gcc.gnu.org/gcc-9/porting_to.html#ompdatasharing for info on the following
// version-specific gcc macro.
#if defined(__GNUC__) && !defined(__clang__) && __GNUC__ < 9
  // gcc-8.X complains about normal syntax. ATS-2 requires different code from all the other
  // machines - these systems must handle const strings differently.
#ifdef draco_isATS2
  // rzansel/sierra
#pragma omp parallel default(none) shared(hostname, std::cout)
#else
  // e.g. ccs-net with gcc-8.3.1 (RHEL 8) / darwin p9 with gcc-8.3.0
#pragma omp parallel default(none) shared(std::cout)
#endif

#else // gcc version >= 9
#pragma omp parallel default(none) shared(num_cpus, hostname, rank, std::cout)
#endif
  {
    int thread = omp_get_thread_num();
    std::string cpuset = rtt_c4::cpuset_to_string(num_cpus);

#pragma omp critical
    {
      std::cout << hostname << " :: Rank " << std::setfill('0') << std::setw(5) << rank
                << ", Thread " << std::setfill('0') << std::setw(3) << thread
                << ", core affinity = " << cpuset << std::endl;
    } // end omp critical
  }   // end omp parallel

  rtt_c4::finalize();
  return (0);
}

//------------------------------------------------------------------------------------------------//
// End c4/bin/xthi.cc
//------------------------------------------------------------------------------------------------//
