//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   c4/test/tstOMP_API_off.cc
 * \author Timothy Kelley
 * \date   Tue Jun  9 15:03:08 2020
 * \brief  Demonstrate basic OMP API.
 * \note   Copyright (C) 2020-2022 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"
#include "ds++/Soft_Equivalence.hh"
#include "ds++/config.h" // OPENMP_FOUND will be set here
#include <functional>    // std::function

#ifdef OPENMP_FOUND
#define ORIG_OPENMP_FOUND OPENMP_FOUND
#undef OPENMP_FOUND
#endif

#include "c4/c4_omp.h"

using namespace rtt_c4;
using rtt_dsxx::UnitTest;

/* If OpenMP is available, we will pretend that it's not, but still use the
 * OpenMP runtime to check that the funtions do not change OpenMP state.
 * If OpenMP truly is not available, we will still check that the bypass
 * versions behave as expected.
 */

#ifdef ORIG_OPENMP_FOUND
#include <omp.h>
#endif

void check_set_get(UnitTest &ut) {
  int const init_n = get_omp_num_threads();
  FAIL_IF_NOT(1 == init_n);

#ifdef ORIG_OPENMP_FOUND
  int const true_init_n = omp_get_max_threads();
#endif

  set_omp_num_threads(51);
  int const new_n = get_omp_num_threads();
  FAIL_IF_NOT(1 == new_n);
  int const new_max = get_omp_max_threads();
  FAIL_IF_NOT(1 == new_max);
  int const thread_num = get_omp_thread_num();
  FAIL_IF_NOT(0 == thread_num);

#ifdef ORIG_OPENMP_FOUND
  // Even if OpenMP is available, the preceding should have changed nothing
  int const true_new_n = omp_get_max_threads();
  FAIL_IF_NOT(true_new_n == true_init_n);
  // If we did make a mess, clean it up
  if (true_new_n != true_init_n) {
    omp_set_num_threads(true_init_n);
  }
#endif

  return;
}

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

using t_func = std::function<void(UnitTest &)>;

void run_a_test(UnitTest &u, t_func f, std::string const &msg) {
  f(u);
  if (u.numFails == 0) {
    u.passes(msg);
  }
  return;
}

//------------------------------------------------------------------------------------------------//
int main(int argc, char **argv) {
  rtt_dsxx::ScalarUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    run_a_test(ut, check_set_get, "Check getting/setting OpenMP number threads, API unavailable");
  } // try--catches in the epilog:
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstOMP_API_off.cc
//------------------------------------------------------------------------------------------------//
