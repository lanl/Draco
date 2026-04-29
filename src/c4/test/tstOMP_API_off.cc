//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   c4/test/tstOMP_API_off.cc
 * \author Timothy Kelley
 * \brief  Demonstrate basic OMP API.
 * \note   Copyright (C) 2020-2026 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "dsxx/Release.hh"
#include "dsxx/ScalarUnitTest.hh"
#include "dsxx/UnitTest.hh"
#include "dsxx/config.h" // OPENMP_FOUND will be set here
#include <string>

#ifdef OPENMP_FOUND
#define ORIG_OPENMP_FOUND OPENMP_FOUND
#undef OPENMP_FOUND
#endif

#include "c4/c4_omp.h"

namespace {

using namespace rtt_c4;
using rtt_dsxx::UnitTest;

// If OpenMP is available, we will pretend that it's not, but still use the OpenMP runtime to check
// that the functions do not change OpenMP state.  If OpenMP truly is not available, we will still
// check that the bypass versions behave as expected.

#ifdef ORIG_OPENMP_FOUND
#include <omp.h>
#endif

void check_set_get(UnitTest &ut) {
  int const init_n = get_omp_num_threads();
  FAIL_IF_NOT(1 == init_n);

#ifdef ORIG_OPENMP_FOUND
  int const true_init_n = omp_get_max_threads();
  constexpr auto kind{omp_sched_guided};
  constexpr int chunk_size{1337};
  omp_set_schedule(kind, chunk_size);
#endif

  set_omp_num_threads(51);
  int const new_n = get_omp_num_threads();
  FAIL_IF_NOT(1 == new_n);
  int const new_max = get_omp_max_threads();
  FAIL_IF_NOT(1 == new_max);
  int const thread_num = get_omp_thread_num();
  FAIL_IF_NOT(0 == thread_num);

  // Static sched value from OMP documentation (we can't use omp_sched_t).
  constexpr auto new_kind{0x1};
  constexpr int new_chunk_size{42};
  set_omp_schedule(new_kind, new_chunk_size);

#ifdef ORIG_OPENMP_FOUND
  // Even if OpenMP is available, the preceding should have changed nothing
  int const true_new_n = omp_get_max_threads();
  FAIL_IF_NOT(true_new_n == true_init_n);
  // If we did make a mess, clean it up
  if (true_new_n != true_init_n) {
    omp_set_num_threads(true_init_n);
  }
  omp_sched_t check_kind;
  int check_chunk;
  omp_get_schedule(&check_kind, &check_chunk);
  // Check that sched. params weren't altered:
  FAIL_IF_NOT(kind == check_kind);
  FAIL_IF_NOT(chunk_size == check_chunk);
  // Clean up:
  set_omp_schedule(omp_sched_auto);
#endif

  return;
}

} // namespace

//------------------------------------------------------------------------------------------------//
int main(int argc, char **argv) {
  using namespace rtt_dsxx;
  ScalarUnitTest ut(argc, argv, release);
  try {
    run_a_test(ut, check_set_get, "Check getting/setting OpenMP number threads, API unavailable");
  } // try--catches in the epilog:
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstOMP_API_off.cc
//------------------------------------------------------------------------------------------------//
