//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   c4/test/tstOMP.cc
 * \author Kelly Thompson
 * \date   Tue Jun  6 15:03:08 2006
 * \brief  Demonstrate basic OMP threads under MPI.
 * \note   Copyright (C) 2011-2021 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "c4/ParallelUnitTest.hh"
#include "c4/Timer.hh"
#include "c4/c4_omp.h"
#include "c4/gatherv.hh"
#include "ds++/Release.hh"
#include "ds++/Soft_Equivalence.hh"
#include <array>
#include <complex>
#include <numeric>

using namespace rtt_c4;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

//------------------------------------------------------------------------------------------------//
bool topology_report() {
  size_t const mpi_ranks = rtt_c4::nodes();
  size_t const my_mpi_rank = rtt_c4::node();

  // Store proc name on local proc
  std::string my_pname = rtt_c4::get_processor_name();
  Remember(size_t namelen = my_pname.size(););

  // Create a container on IO proc to hold names of all nodes.
  std::vector<std::string> procnames(mpi_ranks);

  // Gather names into pnames on IO proc.
  rtt_c4::indeterminate_gatherv(my_pname, procnames);

  // Is there only 1 MPI rank per machine node?
  int one_mpi_rank_per_node(0);

  // Look at the data found on the IO proc.
  if (my_mpi_rank == 0) {
    Check(procnames[my_mpi_rank].size() == namelen);

    // Count unique processors
    std::vector<std::string> unique_processor_names;
    for (size_t i = 0; i < mpi_ranks; ++i) {
      bool found(false);
      for (const auto &unique_processor_name : unique_processor_names)
        if (procnames[i] == unique_processor_name)
          found = true;
      if (!found)
        unique_processor_names.push_back(procnames[i]);
    }

    // Print a report
    std::cout << "\nWe are using " << mpi_ranks << " mpi rank(s) on "
              << unique_processor_names.size() << " unique nodes.";

    for (size_t i = 0; i < mpi_ranks; ++i)
      std::cout << "\n  - MPI rank " << i << " is on " << procnames[i];
    std::cout << "\n" << std::endl;

    if (mpi_ranks == unique_processor_names.size())
      one_mpi_rank_per_node = 1;
  }

  rtt_c4::broadcast(&one_mpi_rank_per_node, 1, 0);

  // return 't' if 1 MPI rank per machine node.
  return (one_mpi_rank_per_node == 1);
}

//------------------------------------------------------------------------------------------------//
void topo_report(rtt_dsxx::UnitTest &ut, bool &one_mpi_rank_per_node) {
  // Determine if MPI ranks are on unique machine nodes:
  //
  // If there are multiple MPI ranks per machine node, then don't use OMP because OMP can't restrict
  // its threads to running only on an MPI rank's cores.  The OMP threads will be distributed over
  // the whole machine node.  For example, we might choose to use 4 MPI ranks on a machine node with
  // 16 cores.  Ideally, we could allow each MPI rank to use 4 OMP threads for a maximum of 4x4=16
  // OMP threads on the 16 core node.  However, because OMP doesn't know about the MPI ranks sharing
  // the 16 cores, the even distribution of OMP threads is not guaranteed.
  //
  // So - if we have more than one MPI rank per machine node, then turn off OMP threads.
  one_mpi_rank_per_node = topology_report();

  std::string procname = rtt_c4::get_processor_name();

#ifdef OPENMP_FOUND

  // Turn on the dynamic thread adjustment capability.
  omp_set_dynamic(1);
  int num_dynamic_threads = omp_get_dynamic();

  int tid(-1);
  int nthreads(-1);
  int maxthreads(-1);

  maxthreads = omp_get_max_threads();
  // This is just a unit test. Limit the parallelism.
  if (maxthreads > 16)
    omp_set_num_threads(16);

#pragma omp parallel private(tid)
  {
    nthreads = omp_get_num_threads();
    tid = omp_get_thread_num();

    if (tid == 0) {
      std::cout << "Using OMP threads."
                << "\n   MPI node       : " << node() << "\n   MPI max nodes  : " << nodes()
                << "\n   OMP thread     : " << tid << "\n   OMP num threads: " << nthreads
                << "\n   OMP max threads: " << maxthreads << "\n   procname(IO)   : " << procname
                << "\n   Dynamic theads : "
                << (num_dynamic_threads == 0 ? std::string("OFF") : std::string("ON")) << "\n"
                << std::endl;
    }
    FAIL_IF(tid < 0 || tid >= nthreads);
  }
#else
  { // not OMP
    std::cout << "OMP thread use is disabled."
              << "\n   MPI node       : " << node() << "\n   MPI max nodes  : " << nodes()
              << "\n   procname(IO)   : " << procname << "\n"
              << std::endl;
    PASSMSG("OMP is disabled.  No checks made.");
  }
#endif

  if (ut.numFails == 0)
    PASSMSG("topology report finished successfully.");
  else
    FAILMSG("topology report failed.");

  return;
}

//------------------------------------------------------------------------------------------------//
void sample_sum(rtt_dsxx::UnitTest &ut, bool const omrpn) {
  if (rtt_c4::node() == 0)
    std::cout << "Begin test sample_sum()...\n" << std::endl;

  // Generate data and benchmark values:
  int N(10000000);
  std::vector<double> foo(N, 0.0);
  std::vector<double> result(N, 0.0);
  std::vector<double> bar(N, 99.0);

  Timer t1_serial_build;
  t1_serial_build.start();

  for (int i = 0; i < N; ++i) {
    foo[i] = 99.00 + i;
    bar[i] = 0.99 * i;
    result[i] = std::sqrt(foo[i] + bar[i]) + 1.0;
  }
  t1_serial_build.stop();

  Timer t2_serial_accumulate;
  t2_serial_accumulate.start();

  double sum = std::accumulate(foo.begin(), foo.end(), 0.0);

  t2_serial_accumulate.stop();

  if (node() == 0)
    std::cout << "benchmark: sum(foo) = " << sum << std::endl;

#ifdef OPENMP_FOUND
  {
    // This is just a unit test. Limit the parallelism.
    int maxthreads(-1);
    maxthreads = omp_get_max_threads();
    if (maxthreads > 16)
      omp_set_num_threads(16);

    // More than 1 MPI rank per node --> turn off OMP.
    if (!omrpn)
      omp_set_num_threads(1);

    // Generate omp_result
    std::vector<double> omp_result(N, 0.0);
    double omp_sum(0.0);

    Timer t1_omp_build;
    t1_omp_build.start();

    int nthreads(-1);
#pragma omp parallel
    {
      if (node() == 0 && omp_get_thread_num() == 0) {
        nthreads = omp_get_num_threads();
        std::cout << "\nNow computing sum using " << nthreads << " OMP threads." << std::endl;
      }
    }

#pragma omp parallel for shared(foo, bar)
    for (int i = 0; i < N; ++i) {
      foo[i] = 99.00 + i;
      bar[i] = 0.99 * i;
      result[i] = std::sqrt(foo[i] + bar[i]) + 1.0;
    }
    t1_omp_build.stop();

    // Accumulate via OMP

    Timer t2_omp_accumulate;
    t2_omp_accumulate.start();

// clang-format adds spaces around this colon.
// clang-format off
#pragma omp parallel for reduction(+: omp_sum)
    // clang-format on
    for (int i = 0; i < N; ++i)
      omp_sum += foo[i];

    t2_omp_accumulate.stop();

    // Sanity check
    if (rtt_dsxx::soft_equiv(sum, omp_sum))
      PASSMSG("OpenMP sum matches std::accumulate() value!");
    else
      FAILMSG("OpenMP sum differs!");

    if (node() == 0) {
      std::cout.precision(6);
      std::cout.setf(std::ios::fixed, std::ios::floatfield);
      std::cout << "Timers:"
                << "\n\t             \tSerial Time \tOMP Time"
                << "\n\tbuild      = \t" << t1_serial_build.wall_clock() << "\t"
                << t1_omp_build.wall_clock() << "\n\taccumulate = \t"
                << t2_serial_accumulate.wall_clock() << "\t" << t2_omp_accumulate.wall_clock()
                << std::endl;
    }

    // [2015-11-17 KT] The accumulate test no longer provides enough work to offset the overhead of
    // OpenMP, especially for the optimized build.  Turn this test off...

    // if( omrpn && nthreads > 4 )
    // {
    //     if( t2_omp_accumulate.wall_clock()
    //         < t2_serial_accumulate.wall_clock() )
    //         PASSMSG( "OMP accumulate was faster than Serial accumulate.");
    //     else
    //         FAILMSG( "OMP accumulate was slower than Serial accumulate.");
    // }
  }
#else // SCALAR
  PASSMSG("OMP is disabled.  No checks made.");
#endif
  return;
}

//------------------------------------------------------------------------------------------------//
// This is a simple demonstration problem for OMP.  Nothing really to check for PASS/FAIL.
int MandelbrotCalculate(std::complex<double> c, int maxiter) {
  // iterates z = z*z + c until |z| >= 2 or maxiter is reached, returns the number of iterations
  std::complex<double> z = c;
  int n = 0;
  for (; n < maxiter; ++n) {
    if (std::abs(z) >= 2.0)
      break;
    z = z * z + c;
  }
  return n;
}

//------------------------------------------------------------------------------------------------//
void MandelbrotDriver(rtt_dsxx::UnitTest &ut) {
  using namespace std;

  const int width = 78;
  const int height = 44;
  const int num_pixels = width * height;
  const complex<double> center(-0.7, 0.0);
  const complex<double> span(2.7, -(4 / 3.0) * 2.7 * height / width);
  const complex<double> begin = center - span / 2.0;
  const int maxiter = 100000;

  // Use OMP threads
  Timer t;
  ostringstream image1;
  ostringstream image2;
  t.start();

  int nthreads(-1);
#ifdef OPENMP_FOUND

  // This is just a unit test. Limit the parallelism.
  int maxthreads(-1);
  maxthreads = omp_get_max_threads();
  if (maxthreads > 16)
    omp_set_num_threads(16);

#pragma omp parallel
  {
    if (node() == 0 && omp_get_thread_num() == 0) {
      nthreads = omp_get_num_threads();
      cout << "\nNow Generating Mandelbrot image (" << nthreads << " OMP threads)...\n" << endl;
    }
  }

#pragma omp parallel for ordered schedule(dynamic)
  for (int pix = 0; pix < num_pixels; ++pix) {
    const int x = pix % width;
    const int y = pix / width;

    complex<double> c =
        begin + complex<double>(x * span.real() / (width + 1.0), y * span.imag() / (height + 1.0));

    int n = MandelbrotCalculate(c, maxiter);
    if (n == maxiter)
      n = 0;

#pragma omp ordered
    {
      char cc = ' ';
      if (n > 0) {
        array<char, 19> const charset = {'.', ',', 'c', '8', 'M', '@', 'j', 'a', 'w', 'r',
                                         'p', 'o', 'g', 'O', 'Q', 'E', 'P', 'G', 'J'};
        cc = charset[n % (charset.size() - 1)];
      }
      image1 << cc;
      if (x + 1 == width)
        image1 << "|\n"; //puts("|");
    }
  }
#endif // OPENMP_FOUND

  t.stop();
  double const gen_time_omp = t.wall_clock();

  // Repeat for serial case
  if (rtt_c4::node() == 0)
    cout << "\nGenerating Mandelbrot image (Serial)...\n" << endl;

  t.reset();
  t.start();

  for (int pix = 0; pix < num_pixels; ++pix) {
    const int x = pix % width;
    const int y = pix / width;

    complex<double> c =
        begin + complex<double>(x * span.real() / (width + 1.0), y * span.imag() / (height + 1.0));

    int n = MandelbrotCalculate(c, maxiter);
    if (n == maxiter)
      n = 0;

    {
      char cc = ' ';
      if (n > 0) {
        array<char, 19> const charset = {'.', ',', 'c', '8', 'M', '@', 'j', 'a', 'w', 'r',
                                         'p', 'o', 'g', 'O', 'Q', 'E', 'P', 'G', 'J'};
        cc = charset[n % (charset.size() - 1)];
      }
      // putchar(c);
      image2 << cc;
      if (x + 1 == width)
        image2 << "|\n"; //puts("|");
    }
  }
  t.stop();
  double const gen_time_serial = t.wall_clock();

#ifdef OPENMP_FOUND
  if (image1.str() == image2.str()) {
    // cout << image1.str() << endl;
    PASSMSG("Scalar and OMP generated Mandelbrot images match.");
  } else {
    FAILMSG("Scalar and OMP generated Mandelbrot images do not match.");
  }
#endif

  cout << "\nTime to generate Mandelbrot:"
       << "\n   Normal: " << gen_time_serial << " sec." << endl;

  if (nthreads > 4) {
    cout << "   OMP   : " << gen_time_omp << " sec." << endl;
    if (gen_time_omp < gen_time_serial)
      PASSMSG("OMP generation of Mandelbrot image is faster.");
    else
      FAILMSG("OMP generation of Mandelbrot image is slower.");
  }

  return;
}

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  rtt_c4::ParallelUnitTest ut(argc, argv, rtt_dsxx::release);
  try {

#ifndef MSVC
    if (rtt_c4::node() == 0) {
      // Machine load:
      CPUinfo().report_load_average();
      CPUinfo().report_cpu_utilization();
    }
#endif

    // One MPI rank per machine node?
    bool omrpn(false);

    // Unit tests
    topo_report(ut, omrpn);
    sample_sum(ut, omrpn);

    if (rtt_c4::nodes() == 1)
      MandelbrotDriver(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstOMP.cc
//------------------------------------------------------------------------------------------------//
