//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   dsxx/AVX_VectorMath.hh
 * \brief  Vectorized math functions for SIMD operations
 * \note   Copyright (C) 2026 Triad National Security, LLC., All rights reserved.
 *
 * This file provides vectorized implementations of common math functions (exp, log, etc.)
 * using AVX2 intrinsics. These are approximations optimized for performance while maintaining
 * reasonable accuracy for physics calculations.
 *
 * This file was generated all or in part by AI.
 */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_dsxx_AVX_VectorMath_hh
#define rtt_dsxx_AVX_VectorMath_hh

#ifdef __AVX2__
#include <cmath>
#include <immintrin.h>

namespace rtt_dsxx {

//================================================================================================//
/*!
 * \brief Vectorized exponential function using AVX2 (processes 4 doubles)
 *
 * This function was generated (nearly entirely) by Google AI
 *
 * \param[in] x Input vector (4 doubles) - assumes x >= 0 for exp(-x)
 * \return exp(-x) for each element
 */
inline __m256d exp_pd_avx2(__m256d xp) {

  // Clamp to prevent overflow (exp(-20) ~ 2e-9, good enough for transport)
  const __m256d max_arg = _mm256_set1_pd(20.0);
  xp = _mm256_min_pd(xp, max_arg);

  // negate xp
  const __m256d z = _mm256_setzero_pd();
  const __m256d x = _mm256_sub_pd(z, xp);

  // Constants for range reduction
  const __m256d ln2_inv = _mm256_set1_pd(1.4426950408889634074); // 1/ln(2)
  const __m256d ln2_hi = _mm256_set1_pd(0.6931471805599453);     // ln(2) parts for precision
  const __m256d ln2_lo = _mm256_set1_pd(2.319042899806505e-17);

  // Polynomial coefficients for e^f where f belongs to [-0.5, 0.5]
  const __m256d c0 = _mm256_set1_pd(1.0);
  const __m256d c1 = _mm256_set1_pd(1.0);
  const __m256d c2 = _mm256_set1_pd(0.5);
  const __m256d c3 = _mm256_set1_pd(0.16666666666666666);
  const __m256d c4 = _mm256_set1_pd(0.041666666666666664);
  const __m256d c5 = _mm256_set1_pd(0.008333333333333333);
  const __m256d c6 = _mm256_set1_pd(0.001388888888888889);

  // 1. Calculate k = round(x / ln(2))
  const __m256d k_m =
      _mm256_round_pd(_mm256_mul_pd(x, ln2_inv), _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC);

  // 2. Subtract k * ln(2) to find the fractional part: f = x - k * ln(2)
  __m256d f = _mm256_fnmadd_pd(k_m, ln2_hi, x);
  f = _mm256_fnmadd_pd(k_m, ln2_lo, f);

  // 3. Evaluate polynomial using Horner's Scheme and FMA3
  __m256d poly = c6;
  poly = _mm256_fmadd_pd(poly, f, c5);
  poly = _mm256_fmadd_pd(poly, f, c4);
  poly = _mm256_fmadd_pd(poly, f, c3);
  poly = _mm256_fmadd_pd(poly, f, c2);
  poly = _mm256_fmadd_pd(poly, f, c1);
  poly = _mm256_fmadd_pd(poly, f, c0);

  // 4. Compute 2^k via bit-shifting into the exponent field
  // Convert double precision format to 64-bit integers
  const __m128i k_int = _mm256_cvtpd_epi32(k_m); // Cast down to 32-bit int elements

  // Shift integers to position them into the double-precision exponent layout
  // Double bias is 1023
  const __m128i bias = _mm_set1_epi32(1023);
  const __m128i k_biased = _mm_add_epi32(k_int, bias);

  // Unpack 32-bit ints to 64-bit ints to align with double lanes
  const __m256i k_64 = _mm256_cvtepi32_epi64(k_biased);
  const __m256i twok_bits = _mm256_slli_epi64(k_64, 52);
  const __m256d twok = _mm256_castsi256_pd(twok_bits);

  // 5. Combine: result = poly * 2^k
  return _mm256_mul_pd(poly, twok);
}

//================================================================================================//
/*!
 * \brief Horizontal maximum of a vector (returns the maximum of 4 values)
 *
 * \param[in] v Input vector (4 doubles)
 * \return Maximum value across all 4 elements
 */
inline double horizontal_max_pd(__m256d v) {
  // Swap high and low 128-bit lanes and take max
  const __m256d v_perm = _mm256_permute2f128_pd(v, v, 0x01);
  __m256d v_max = _mm256_max_pd(v, v_perm);

  // Swap adjacent pairs and take max
  const __m256d v_shuffle = _mm256_shuffle_pd(v_max, v_max, 0x05);
  v_max = _mm256_max_pd(v_max, v_shuffle);

  // Extract the maximum value
  return _mm256_cvtsd_f64(v_max);
}

} // end namespace rtt_dsxx

#endif // __AVX2__

#endif // rtt_dsxx_AVX_VectorMath_hh

//------------------------------------------------------------------------------------------------//
// end of dsxx/AVX_VectorMath.hh
//------------------------------------------------------------------------------------------------//
