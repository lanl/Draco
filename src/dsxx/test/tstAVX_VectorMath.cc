//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   dsxx/test/tstAVX_VectorMath.cc
 * \brief  Test functions defined in dsxx/AVX_VectorMath.hh.
 * \note   Copyright (C) 2026 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "dsxx/Release.hh"
#include "dsxx/ScalarUnitTest.hh"
#include "dsxx/UnitTest.hh"

// Include AVX2 support for vectorized streaming
#ifdef __AVX2__
#include "dsxx/AVX_VectorMath.hh"
#include "dsxx/Soft_Equivalence.hh"
#include <algorithm>
#include <array>
#include <cmath>
#include <immintrin.h>
#endif

namespace {

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//
// test implementation of exp
void test_exp_pd_avx2(rtt_dsxx::UnitTest &ut) {

#ifdef __AVX2__
  // test exp_pd_avx2 for 0 through 19
  std::array<double, 4> x = {0.0, 1.0, 2.0, 3.0};
  for (int i = 0; i < 20; i += 4) {
    const __m256d x_vec = _mm256_loadu_pd(&x[0]);
    const __m256d e_vec = rtt_dsxx::exp_pd_avx2(x_vec);
    for (int j = 0; j < 4; ++j)
      FAIL_IF_NOT(rtt_dsxx::soft_equiv(e_vec[j], std::exp(-x[j]), 1e-6));

    // increment x
    std::ranges::transform(x, x.begin(), [](double xv) noexcept { return xv + 4.0; });
  }

  if (ut.numFails == 0)
    PASSMSG("exp_pd_avx2 tests okay.");
  return;
#else
  if (ut.numFails == 0)
    PASSMSG("exp_pd_avx2 (AVX2) not supported for this build -- skipping test.");
  return;
#endif
}

// test avx vector max function
void test_horizontal_max_pd(rtt_dsxx::UnitTest &ut) {

#ifdef __AVX2__
  std::array<double, 4> x = {3.0, 0.0, 1.0, 2.0};
  for (int i = 0; i < 3; ++i) {
    const __m256d x_vec = _mm256_loadu_pd(&x[0]);
    const double mval = rtt_dsxx::horizontal_max_pd(x_vec);
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(mval, x[i]));

    // shift x
    std::ranges::rotate(x.begin(), x.begin() + 3, x.end());
  }

  if (ut.numFails == 0)
    PASSMSG("horizontal_max_pd tests okay.");
  return;
#else
  if (ut.numFails == 0)
    PASSMSG("horizontal_max_pd (AVX2) not supported for this build -- skipping test.");
  return;
#endif
}

} // namespace

int main(int argc, char *argv[]) {
  rtt_dsxx::ScalarUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    test_exp_pd_avx2(ut);
    test_horizontal_max_pd(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstAVX_VectorMath.cc
//------------------------------------------------------------------------------------------------//
