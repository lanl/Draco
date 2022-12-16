//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   ds++/test/tstPrefetch.cc
 * \brief  Demonstrate/Test the prefetch function.
 * \note   Copyright (C) 2018-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "ds++/Prefetch.hh"
#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"
#include "ds++/Soft_Equivalence.hh"
#include <cmath>
#include <iostream>

using namespace std;
using namespace rtt_dsxx;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

void runtest(rtt_dsxx::UnitTest &ut) {

  cout << "Begin tstPrefetch::runtest() checks...\n";

  // Loop without prefetch
  unsigned const R = 300;

  // Reduce runtime for Debug builds.
#ifdef DEBUG
  unsigned const N = 10000;
#else
  unsigned const N = 1000000;
#endif
  vector<double> huge(N);
  double sum = 0.0;
  for (unsigned k = 0; k < R; ++k) {
    for (unsigned i = 0; i < N; i += CACHE_LINE_DOUBLE) {
      for (unsigned j = 0; j < CACHE_LINE_DOUBLE; ++j) {
        huge[j] = 0.12 * j;
      }
    }
    for (unsigned i = 0; i < N; i += CACHE_LINE_DOUBLE) {
      for (size_t j = 0; j < static_cast<size_t>(CACHE_LINE_DOUBLE); j += 2) {
        sum += sqrt(huge[j] * huge[j] + huge[j + 1] * huge[j + 1]);
      }
    }
  }

  // Save and reset
  double const sum_noprefetch(sum);
  sum = 0.0;

  // Loop with prefetch
  for (unsigned k = 0; k < R; ++k) {
    for (unsigned i = 0; i < N; i += CACHE_LINE_DOUBLE) {
      prefetch_cache_line(&huge[i + CACHE_LINE_DOUBLE], 1, 0);
      for (unsigned j = 0; j < CACHE_LINE_DOUBLE; ++j) {
        huge[j] = 0.12 * j;
      }
    }
    for (unsigned i = 0; i < N; i += CACHE_LINE_DOUBLE) {
      prefetch_cache_line(&huge[i + CACHE_LINE_DOUBLE], 0, 1);
      for (size_t j = 0; j < static_cast<size_t>(CACHE_LINE_DOUBLE); j += 2) {
        sum += sqrt(huge[j] * huge[j] + huge[j + 1] * huge[j + 1]);
      }
    }
  }
  FAIL_IF_NOT(soft_equiv(sum_noprefetch, sum));
  PASSMSG("run to completion"); // Nothing crashed
}

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  rtt_dsxx::ScalarUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    runtest(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
//  end of tstPrefetch.cc
//------------------------------------------------------------------------------------------------//
