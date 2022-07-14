//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   c4/test/tstQuoWrapper.cc
 * \author Kelly Thompson
 * \date   Friday, Nov 29, 2019, 19:48 pm
 * \brief  C4 QuoWrapper test.
 * \note   Copyright (C) 2019-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "c4/ParallelUnitTest.hh"
#include "c4/QuoWrapper.hh"
#include "c4/SLURM_Task_Info.hh"
#include "c4/bin/ythi.hh"
#include "ds++/DracoMath.hh"
#include "ds++/Release.hh"
#include <chrono>

//------------------------------------------------------------------------------------------------//
// Helper functions
//------------------------------------------------------------------------------------------------//

// Add barriers and sleep_for in an attempt to keep the output synchronized between ranks and
// threads.
void sync_output() {
  std::cout.flush();
#ifdef C4_MPI
#ifdef HAVE_LIBQUO
  QUO_barrier(rtt_c4::QuoWrapper::fetch_quo());
  rtt_c4::global_barrier();
#endif
#endif
  std::this_thread::sleep_for(std::chrono::seconds(1));
  return;
}

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

void quo_hw_report(rtt_dsxx::UnitTest &ut) {
  using namespace std;

  if (rtt_c4::rank() == 0)
    cout << "\n>>> Starting the quo_hw_report test..." << endl;

    // Only run these checks if MPI is available.
#ifdef C4_MPI

    // Only run these checks if MPI and Libquo are available.
#ifdef HAVE_LIBQUO

  // Generic checks.
  FAIL_IF_NOT(rtt_c4::QuoWrapper::num_nodes() > 0);
  FAIL_IF_NOT(rtt_c4::QuoWrapper::num_cores() > 0);

  // Qup will give the available system resources in the num_hw_threads call and
  // std::thread::hardware_concurency will give the total system resources (not respecting what has
  // been reserved for the OS). Use std::thread as an upper bound check
  FAIL_IF_NOT(rtt_c4::QuoWrapper::num_hw_threads() <= thread::hardware_concurrency());

  FAIL_IF_NOT(rtt_c4::QuoWrapper::num_sockets_per_node() > 0);
  FAIL_IF_NOT(rtt_c4::QuoWrapper::num_numanodes_per_node() > 0);
  FAIL_IF_NOT(rtt_c4::QuoWrapper::num_mpi_ranks_per_node() > 0);
  FAIL_IF_NOT(rtt_c4::QuoWrapper::num_mpi_ranks_per_node() <= rtt_c4::nranks());
  FAIL_IF_NOT(rtt_c4::QuoWrapper::bindings().length() > 0);

  bool using_hyperthreads =
      (rtt_c4::QuoWrapper::num_cores() != rtt_c4::QuoWrapper::num_hw_threads()) ? true : false;
  if (rtt_c4::rank() == 0) {
    cout << "\nThe current application has access to the following resources:"
         << "\n - Nodes used by this process : " << rtt_c4::QuoWrapper::num_nodes()
         << "\n - Cores per node             : " << rtt_c4::QuoWrapper::num_cores()
         << "\n - Hardware threads per node  : " << rtt_c4::QuoWrapper::num_hw_threads();
    if (using_hyperthreads) {
      cout << " (hyperthreading enabled)";
    }
    cout << "\n - Sockets per node           : " << rtt_c4::QuoWrapper::num_sockets_per_node()
         << "\n - NumaNodes per node         : " << rtt_c4::QuoWrapper::num_numanodes_per_node()
         << "\n - MPI ranks on this node     : " << rtt_c4::QuoWrapper::num_mpi_ranks_per_node()
         << endl
         << endl;
  }
  cout << "MPI Rank " << rtt_c4::rank() << " - Bindings: " << rtt_c4::QuoWrapper::bindings()
       << endl;
  sync_output();

  if (rtt_c4::rank() == 0)
    cout << endl;

  // If we have access to slurm, add checks to compare Libquo values to slurm values.  \c
  // is_job_num_nodes_set will only be true if SLURM is available.
  {
    rtt_c4::SLURM_Task_Info sti;
    if (sti.is_job_num_nodes_set()) {

      // get_job_num_nodes == max nodes available to this process.
      // num_nodes         == num nodes actually used by this process.
      FAIL_IF_NOT(rtt_c4::QuoWrapper::num_nodes() <= sti.get_job_num_nodes());
      // SLURM_CPUS_ON_NODE is not set consistently across all machines so it is nearly impossible
      // to compare Quo values to SLURM values.  Use inequality checks instead.
      FAIL_IF_NOT(rtt_c4::QuoWrapper::num_hw_threads() >= sti.get_cpus_on_node());
      FAIL_IF_NOT(rtt_c4::QuoWrapper::num_cores() <= sti.get_cpus_on_node());
      FAIL_IF_NOT(rtt_c4::QuoWrapper::num_mpi_ranks_per_node() <= sti.get_cpus_on_node());

      if (rtt_c4::rank() == 0) {
        cout << "\nSLURM values:"
             << "\n - cpus per task (-c) = ";
        if (sti.get_cpus_per_task() < 10000000)
          cout << sti.get_cpus_per_task();
        else
          cout << "(not set, assume 1)";
        cout << "\n - ntasks (-n)        = " << sti.get_ntasks()
             << "\n - nnodes (-N)        = " << sti.get_job_num_nodes()
             << "\n - cpus per node      = " << sti.get_cpus_on_node()
             << "\n - nodelist           = " << sti.get_nodelist() << endl;
      }
    }
  }
  rtt_c4::global_barrier();

  if (rtt_c4::rank() == 0)
    cout << endl;

#else
  PASSMSG("Quo not found. Tests skipped.");
#endif

    // scalar mode ----------------------------------------
#else

#ifdef HAVE_LIBQUO
  // Should this throw or return zero?
  FAIL_IF_NOT(rtt_c4::QuoWrapper::num_cores() == 0);
  FAIL_IF_NOT(rtt_c4::QuoWrapper::num_hw_threads() == 0);
  FAIL_IF_NOT(rtt_c4::QuoWrapper::num_mpi_ranks_per_node() == 0);
  FAIL_IF_NOT(rtt_c4::QuoWrapper::num_sockets_per_node() == 0);
  FAIL_IF_NOT(rtt_c4::QuoWrapper::num_numanodes_per_node() == 0);
  FAIL_IF_NOT(rtt_c4::QuoWrapper::bindings().length() == 0);

#else
  PASSMSG("Quo not found. Tests skipped.");
#endif

#endif

  rtt_c4::global_barrier();
  if (rtt_c4::rank() == 0)
    PASSMSG("Done with tests that querying Quo environment.");
  return;
}

//------------------------------------------------------------------------------------------------//
void quo_new_bindings(rtt_dsxx::UnitTest &ut) {
  using namespace std;

  if (rtt_c4::rank() == 0)
    cout << "\n>>> Starting the quo_new_bindings test...\n" << endl;

#ifdef C4_MPI
#ifdef HAVE_LIBQUO

  // handle to the quo context.
  QUO_t *quo = rtt_c4::QuoWrapper::fetch_quo();

  // Target resource: NUMANODE
  // \sa https://github.com/lanl/libquo/blob/master/src/quo.h
  QUO_obj_type_t const QUO_resource = {QUO_OBJ_NUMANODE};

  // maximum number of processes per resource (e.g.: processes per numanode).  (32 cores / node ) *
  // (1 node / 4 numanodes) = 8 processes/numanode
  uint32_t const nmrpn = rtt_c4::QuoWrapper::num_mpi_ranks_per_node();
  uint32_t const nnnpn = rtt_c4::QuoWrapper::num_numanodes_per_node();
  uint32_t const nproc_per_res = rtt_dsxx::ceil_int_division(nmrpn, nnnpn);

  // number of threads per MPI rank to fill all available cores (resource == numanode)
  //
  // To use all physical resources, total num threads /node = num cores / node.
  //
  // e.g.: (32 cores / node) * ( 1 node / 8 mpi ranks ) = 4 threads / mpi rank.
  uint32_t const nc = rtt_c4::QuoWrapper::num_cores();
  // limit max num_workers to (num_cores/numa -1) to avoid threads talking across a numa
  // boundary. Note: num_workers = threads/mpi-rank - 1;
  uint32_t const num_workers =
      std::min(rtt_dsxx::ceil_int_division(nc, nmrpn), rtt_dsxx::ceil_int_division(nc, nnnpn)) - 1;

  if (rtt_c4::node() == 0) {
    cout << "ncores per node = " << rtt_c4::QuoWrapper::num_cores()
         << "\nmpi ranks       = " << rtt_c4::nranks()
         << "\nnproc_per_res (mpi ranks / numa) = " << nproc_per_res
         << "\nnum_workers (threads / mpi rank) = " << num_workers
         << " (target value, +1 for threads/mpi-rank)"
         << "\n"
         << endl;
  }

  if (rtt_c4::rank() == 0)
    cout << "\n>> Initial bindings and thread affinity...\n" << endl;

  // print a report of mpi-rank and thread bindings
  rtt_c4::report_bindings(num_workers);
  sync_output();

  if (rtt_c4::rank() == 0)
    cout << "\n>> Call Quo_auto_distrib and push new Quo context...\n" << endl;
  sync_output();

  // let quo distribute workers over the sockets, if res_assigned is 1 after this call then I have
  // been chosen. (This isn't needed if MPI bindings are not changing.)
  //
  // int res_assigned(0);
  // QUO_auto_distrib(quo, QUO_resource, nproc_per_res, &res_assigned);
  // if (res_assigned != 0) {

  // Push the new thread bindings onto the stack (activate)
  QUO_bind_push(quo, QUO_BIND_PUSH_OBJ, QUO_resource, 0);

  rtt_c4::report_bindings(num_workers);

  QUO_barrier(quo);
  QUO_bind_pop(quo);
  sync_output();

  //}

  if (rtt_c4::rank() == 0)
    cout << "\n>> Return to original bindings...\n" << endl;
  sync_output();

  rtt_c4::report_bindings(num_workers);
  sync_output();

#endif
#endif // C4_MPI

  // no C4_Serial tests...
  // no tests if QUO not available.

  rtt_c4::global_barrier();
  if (rtt_c4::rank() == 0)
    PASSMSG("Done with tests that modify MPI and thread layout..");
  return;
}

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  rtt_c4::ParallelUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    quo_hw_report(ut);
    quo_new_bindings(ut);
    std::cout << "==> done with tests. cleaning up..." << std::endl;
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstQuoWrapper.cc
//------------------------------------------------------------------------------------------------//
