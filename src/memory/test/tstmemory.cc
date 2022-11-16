//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   memory/test/tstmemory.cc
 * \author Kent G. Budge, Kelly G. Thompson
 * \brief  memory test.
 * \note   Copyright (C) 2013-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"
#include "memory/memory.hh"
#include <sstream>

using namespace std;
using namespace rtt_memory;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

size_t get_really_big_size_t();

void tst_memory(rtt_dsxx::UnitTest &ut) {
  if (total_allocation() == 0) {
    PASSMSG("correct initial total allocation");
  } else {
    FAILMSG("NOT correct initial total allocation");
  }
  if (peak_allocation() == 0) {
    PASSMSG("correct initial peak allocation");
  } else {
    FAILMSG("NOT correct initial peak allocation");
  }
  if (largest_allocation() == 0) {
    PASSMSG("correct initial largest allocation");
  } else {
    FAILMSG("NOT correct initial largest allocation");
  }

  set_memory_checking(true);

  auto *array = new double[20];
  auto *array2 = new double[30];

#if DRACO_DIAGNOSTICS & 2
  if (total_allocation() == 50 * sizeof(double)) {
    PASSMSG("correct total allocation");
  } else {
    FAILMSG("NOT correct total allocation");
  }
  if (peak_allocation() >= 50 * sizeof(double)) {
    PASSMSG("correct peak allocation");
  } else {
    FAILMSG("NOT correct peak allocation");
  }
  if (largest_allocation() >= 30 * sizeof(double)) {
    PASSMSG("correct largest allocation");
  } else {
    FAILMSG("NOT correct largest allocation");
  }
#else
  PASSMSG("memory diagnostics not checked for this build");
#endif
  report_leaks(cerr);

  delete[] array;
  delete[] array2;

#if DRACO_DIAGNOSTICS & 2
  if (total_allocation() == 0) {
    PASSMSG("correct total allocation");
  } else {
    FAILMSG("NOT correct total allocation");
  }
  if (peak_allocation() >= 50 * sizeof(double)) {
    PASSMSG("correct peak allocation");
  } else {
    FAILMSG("NOT correct peak allocation");
  }
  if (largest_allocation() >= 30 * sizeof(double)) {
    PASSMSG("correct largest allocation");
  } else {
    FAILMSG("NOT correct largest allocation");
  }
#endif

  // Just to try to exercise the sized delete version.
  int *scalar = new int;

#if DRACO_DIAGNOSTICS & 2
  if (total_allocation() == sizeof(int)) {
    PASSMSG("correct total allocation");
  } else {
    FAILMSG("NOT correct total allocation");
  }
  if (peak_allocation() >= 50 * sizeof(double) + sizeof(int)) {
    PASSMSG("correct peak allocation");
  } else {
    FAILMSG("NOT correct peak allocation");
  }
  if (largest_allocation() >= 30 * sizeof(double)) {
    PASSMSG("correct largest allocation");
  } else {
    FAILMSG("NOT correct largest allocation");
  }
#endif

  delete scalar;

#if DRACO_DIAGNOSTICS & 2
  if (total_allocation() == 0) {
    PASSMSG("correct total allocation");
  } else {
    FAILMSG("NOT correct total allocation");
  }
  if (peak_allocation() >= 50 * sizeof(double) + sizeof(int)) {
    PASSMSG("correct peak allocation");
  } else {
    FAILMSG("NOT correct peak allocation");
  }
  if (largest_allocation() >= 30 * sizeof(double)) {
    PASSMSG("correct largest allocation");
  } else {
    FAILMSG("NOT correct largest allocation");
  }
  report_leaks(cerr);
#endif
}

//------------------------------------------------------------------------------------------------//
void tst_bad_alloc(rtt_dsxx::UnitTest &ut) {
  size_t const num_fails_start(ut.numFails);

  std::cout << "\nTesting special handler (stack trace) for std::bad_alloc...\n" << std::endl;

  // Set a specialized memory handler.
  std::set_new_handler(rtt_memory::out_of_memory_handler);

  try {
    // trigger a std::bad_alloc exception
    std::cout << "Attempt to allocate some memory." << std::endl;

    char *pBigArray = new char[get_really_big_size_t()];

    // should never get here.
    {
      pBigArray[0] = 'a';
      std::cout << "pBigArray[0] = " << pBigArray[0]
                << "\ntotal_allocation = " << total_allocation() << std::endl;
    }

    delete[] pBigArray;
  } catch (std::bad_alloc & /*error*/) {
    std::cout << "Successfully caught an expected std::bad_alloc exception." << std::endl;
  } catch (...) {
    FAILMSG("Failed to catch a bad_alloc.");
  }

  if (ut.numFails > num_fails_start)
    FAILMSG("Test failures in tst_bad_alloc.");
  else
    PASSMSG("tst_bad_alloc completed successfully.");

  return;
}

//------------------------------------------------------------------------------------------------//
void tst_report(rtt_dsxx::UnitTest & /*unused*/) {
  // No good way to test automatically, alas. All we check is that the calls go in and don't crash
  // anything.
  {
    vector<unsigned char> d1(1000);
    rtt_memory::set_report_threshold(1500);
    d1.resize(2000);
    d1.clear();
    rtt_memory::set_report_threshold(100);
    d1.resize(3000);
  }
  rtt_memory::set_report_threshold();
}

//------------------------------------------------------------------------------------------------//
// Some compilers are clever enough to figure out that if you pass
// std::numeric_limits<size_t>::max() to the new operator, you will always blow away member, and so
// they will refuse to compile the code. We have to use a bit of indirection to get such compilers
// to swallow the huge allocation meant to deliberately blow away memory.
size_t get_really_big_size_t() { return numeric_limits<std::ptrdiff_t>::max() / sizeof(size_t); }

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  rtt_dsxx::ScalarUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    tst_memory(ut);
    tst_report(ut);
    tst_bad_alloc(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstmemory.cc
//------------------------------------------------------------------------------------------------//
