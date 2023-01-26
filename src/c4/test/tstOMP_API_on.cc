//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   c4/test/tstOMP_API_on.cc
 * \author Timothy Kelley
 * \date   Tue Jun  9 15:03:08 2020
 * \brief  Demonstrate basic OMP API with OpenMP present
 * \note   Copyright (C) 2020-2023 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "c4/c4_omp.h"
#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"
#include "ds++/Soft_Equivalence.hh"
#include "ds++/config.h" // OPENMP_FOUND will be set here
#include <functional>    // std::function

#ifdef OPENMP_FOUND

using namespace rtt_c4;
using rtt_dsxx::UnitTest;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//
void check_set_get(UnitTest &ut) {
  int const init_n = get_omp_num_threads();
  set_omp_num_threads(39);
  // use the direct OMP interface to check number of threads was set correctly
  int const new_max{omp_get_max_threads()};
  FAIL_IF_NOT(39 == new_max);
  int const new_max_us(get_omp_max_threads());
  FAIL_IF_NOT(39 == new_max_us);
  // now reset to the previous number of threads
  set_omp_num_threads(init_n);
  int const final_n{get_omp_num_threads()};
  FAIL_IF_NOT(final_n == init_n);
  int const thread_num{get_omp_thread_num()};
  FAIL_IF(thread_num >= final_n || thread_num < 0);
  // set the schedule to dynamic (see documentation):
  constexpr omp_sched_t kind{omp_sched_dynamic};
  constexpr int chunk_size{42};
  set_omp_schedule(kind, chunk_size);
  // create return vals:
  omp_sched_t kind_check;
  int chunk_check;
  // check that get matches what was set:
  get_omp_schedule(kind_check, chunk_check);
  FAIL_IF_NOT(kind == kind_check);
  FAIL_IF_NOT(chunk_size == chunk_check);
  return;
} // check_set_get

using t_func = std::function<void(UnitTest &)>;

void run_a_test(UnitTest &u, t_func const &f, std::string const &msg) {
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
    run_a_test(ut, check_set_get, "Check getting/setting OpenMP number threads, API available");
  } // try--catches in the epilog:
  UT_EPILOG(ut);
}

#else // OPENMP_FOUND

/* No OpenMP case handled in tstOMP_API_off */
int main(int argc, char **argv) {
  rtt_dsxx::ScalarUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    PASSMSG("OpenMP disabled. No testing will be done.");
  }
  UT_EPILOG(ut);
}

#endif // OPENMP_FOUND

//------------------------------------------------------------------------------------------------//
// end of tstOMP_API_on.cc
//------------------------------------------------------------------------------------------------//
