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

#if defined(__GNUC__) && __GNUC__ == 8 && __GNUC_MINOR__ == 3 && __GNUC_PATCHLEVEL__ == 1
  // gcc-8.3.1 complains about normal syntax. power9 and x86_64 have different requirements.
#ifdef draco_isPPC
  // rzansel/sierra
#pragma omp parallel default(none) shared(hostname, std::cout)
#else
  // ccs-net with gcc-8.3.1 (RHEL 8)
#pragma omp parallel default(none) shared(std::cout)
#endif

#elif defined(__GNUC__) && __GNUC__ == 8 && __GNUC_MINOR__ == 5 && __GNUC_PATCHLEVEL__ == 0
  // ccs-net with gcc-8.5.0 (RHEL 8)
#pragma omp parallel default(none) shared(std::cout)

#else
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
