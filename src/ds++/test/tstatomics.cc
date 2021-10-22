//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   ds++/test/tstatomics.cc
 * \author Tim Kelley
 * \date   Thursday, Sept. 6, 2018, 10:51 am
 * \note   Copyright (C) 2018-2021 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"
#include "ds++/Soft_Equivalence.hh"
#include "ds++/atomics.hh"
#include <functional> // std::function
#include <thread>

using rtt_dsxx::UnitTest;

//------------------------------------------------------------------------------------------------//
/* Hammer an atomic from each thread. Each iteration, the thread adds (tid * iteration) to the
 * counter. The atomic ensures that everyone sees a consistent view of the counter: no thread
 * overwrites the contribution from any other thread.
 */
void thread_action(std::atomic<double> &d, size_t N, size_t tid) {
  auto const did = static_cast<double>(tid);
  double d_i = 1;
  for (size_t i = 0; i < N; ++i) {
    double addend = did * d_i;
    rtt_dsxx::fetch_add(d, addend);
    d_i += 1.0;
  }
  return;
} // thread_action

//------------------------------------------------------------------------------------------------//
/* Test fetch_add using an atomic. Expect to get the correct sum every time.*/
void fetch_add_atomic_core(UnitTest &ut, size_t const n_threads, size_t const n_iterations) {
  std::atomic<double> a_d(0.0);

  // launch a number of threads
  std::vector<std::thread> threads(n_threads);
  size_t tid(0);
  for (auto &t : threads) {
    t = std::thread(thread_action, std::ref(a_d), n_iterations, tid);
    tid++;
  }
  // join
  for (auto &t : threads) {
    t.join();
  }

  // Compute expected value the easy way
  double result = a_d.load();
  double sum = 0.0;
  for (size_t i = 0; i < n_iterations; ++i) {
    sum += static_cast<double>(i + 1);
  }
  double tsum = 0.0;
  for (size_t t = 0; t < n_threads; ++t) {
    tsum += static_cast<double>(t);
  }
  double expected = sum * tsum;

  // check and report
  bool const passed = rtt_dsxx::soft_equiv(result, expected);
  if (!passed) {
    printf("%s:%i tsum = %.0f, isum = %.0f, result = %.0f\n", __FUNCTION__, __LINE__, tsum, sum,
           result);
  }
  FAIL_IF_NOT(passed);
  return;
} // fetch_add_atomic_core

void test_fetch_add_atomic(UnitTest &ut) {
  size_t const n_threads(19);
  size_t const n_iterations(10001);
  fetch_add_atomic_core(ut, n_threads, n_iterations);
  return;
} // test_fetch_add_atomic

void test_fetch_add_atomic_1e6(UnitTest &ut) {
  size_t const n_threads(19);
  size_t const n_iterations(1000001);
  fetch_add_atomic_core(ut, n_threads, n_iterations);
  return;
} // test_fetch_add_atomic

// --------------------- non-atomic version --------------------------
// This should give the wrong answer nearly every time on any respectable thread implementation.

//------------------------------------------------------------------------------------------------//
/* Similarly, hammer a POD from each thread. Each iteration, the thread adds (tid * iteration) to
 * the counter. Since the threads are contending, we expect to have a race condition where two
 * threads read the same value from d and one of the thread's write (+=) overwrites the other's.
 */
void thread_action_pod(double &d, size_t N, size_t tid) {
  auto const did = static_cast<double>(tid);
  double d_i = 1;
  for (size_t i = 0; i < N; ++i) {
    double addend = did * d_i;
    d += addend;
    d_i += 1.0;
  }
  return;
} // run_in_a_thread_d

//------------------------------------------------------------------------------------------------//
// same as above, except does not use an atomic
void test_fetch_add_not_atomic(UnitTest & /*ut*/) {
  size_t const n_threads(43);
  size_t const n_iterations(10001);
  double a_d(0.0);

  // launch a number of threads
  std::vector<std::thread> threads(n_threads);
  size_t tid(0);
  for (auto &t : threads) {
    t = std::thread(thread_action_pod, std::ref(a_d), n_iterations, tid);
    tid++;
  }
  // join
  for (auto &t : threads) {
    t.join();
  }

  // calculate expected value
  double result = a_d;
  double sum = 0.0;
  for (size_t i = 0; i < n_iterations; ++i) {
    sum += static_cast<double>(i + 1);
  }
  double tsum = 0.0;
  for (size_t t = 0; t < n_threads; ++t) {
    tsum += static_cast<double>(t);
  }
  double expected = sum * tsum;
  // check and report
  bool const passed = !rtt_dsxx::soft_equiv(result, expected);
  if (!passed) {
    double diff = (expected - result);
    double rel_diff = 100 * diff / expected;
    printf("%s:%i Expected these to differ: tsum = %.0f, isum = %.0f, result = "
           "%.0f, diff = %.0f, rel. "
           "diff = %.2f %% \n",
           __FUNCTION__, __LINE__, tsum, sum, result, diff, rel_diff);
  }
  /* This does not fail on all platforms: on 4 April 2019 failed on appVeyor CI.
   * So, hmmm... */
  // FAIL_IF_NOT(passed);
  return;
} // test_fetch_add_not_atomic

// fetch_sub tests

/* Same as thread_action above, except uses fetch_sub. Total sum is just the negative of the
 * preceding test.
 */
void thread_action_sub(std::atomic<double> &d, size_t N, size_t tid) {
  auto const did = static_cast<double>(tid);
  double d_i = 1;
  for (size_t i = 0; i < N; ++i) {
    double addend = did * d_i;
    rtt_dsxx::fetch_sub(d, addend);
    d_i += 1.0;
  }
  return;
} // thread_action

//------------------------------------------------------------------------------------------------//
/* Test fetch_add using an atomic. Expect to get the correct sum every time.*/
void fetch_sub_atomic_core(UnitTest &ut, size_t const n_threads, size_t const n_iterations) {
  std::atomic<double> a_d(0.0);

  // launch a number of threads
  std::vector<std::thread> threads(n_threads);
  size_t tid(0);
  for (auto &t : threads) {
    t = std::thread(thread_action_sub, std::ref(a_d), n_iterations, tid);
    tid++;
  }
  // join
  for (auto &t : threads) {
    t.join();
  }

  // Compute expected value the easy way
  double result = a_d.load();
  double sum = 0.0;
  for (size_t i = 0; i < n_iterations; ++i) {
    sum += static_cast<double>(i + 1);
  }
  double tsum = 0.0;
  for (size_t t = 0; t < n_threads; ++t) {
    tsum += static_cast<double>(t);
  }
  double expected = -1.0 * sum * tsum;

  // check and report
  bool const passed = rtt_dsxx::soft_equiv(result, expected);
  if (!passed) {
    printf("%s:%i tsum = %.0f, isum = %.0f, result = %.0f, "
           "expected = %.0f\n",
           __FUNCTION__, __LINE__, tsum, sum, result, expected);
  }
  FAIL_IF_NOT(passed);
  return;
} // fetch_add_atomic_core

void test_fetch_sub_atomic(UnitTest &ut) {
  size_t const n_threads(8);
  size_t const n_iterations(10001);
  fetch_sub_atomic_core(ut, n_threads, n_iterations);
  return;
} // test_fetch_add_atomic

void test_fetch_sub_atomic_1e6(UnitTest &ut) {
  size_t const n_threads(8);
  size_t const n_iterations(1000001);
  fetch_sub_atomic_core(ut, n_threads, n_iterations);
  return;
} // test_fetch_add_atomic

using t_func = std::function<void(UnitTest &)>;

//------------------------------------------------------------------------------------------------//
void run_a_test(UnitTest &u, t_func f, std::string const &msg) {
  f(u);
  if (u.numFails == 0) {
    u.passes(msg);
  }
  return;
}

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  rtt_dsxx::ScalarUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    run_a_test(ut, test_fetch_add_atomic, "fetch_add ok.");
    run_a_test(ut, test_fetch_add_atomic_1e6, "fetch_add ok 1e6 iterations.");
    run_a_test(ut, test_fetch_sub_atomic, "fetch_sub ok.");
    run_a_test(ut, test_fetch_sub_atomic_1e6, "fetch_sub ok 1e6 iterations.");
    run_a_test(ut, test_fetch_add_not_atomic, "non-atomic as expected.");
  } // try--catches in the epilog:
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of ds++/test/tstatomics.cc
//------------------------------------------------------------------------------------------------//
