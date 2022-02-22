//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file    ds++/FMA.hh
 * \author  Kelly Thompson
 * \date    Thursday, Feb 09, 2017, 11:22 am
 * \brief   Provide extra control for FMA operations.
 * \note    Copyright (C) 2017-2022 Triad National Security, LLC., All rights reserved.
 *
 * Intel Haswell and later (and also modern AMD cpus) have hardware FMA features. On machines
 * without hardware FMA, the use of \c fma(c,b,a) provides a accurate value than \c a*b+c because
 * roundoff error is accounted for.  However, \c fma can be slow (maybe 3x slower).
 *
 * To improve solution consistency across platforms, we can choose to use a software-based infinite
 * precision fma call.
 *
 * \note Intel's code for detecting FMA availability on hardware.
 * https://software.intel.com/en-us/node/405250?language=es&wapkw=avx2+cpuid
 *
 * \note The Intel compiler doesn't seem to set \c FP_FAST_FMA. If you know that your algorithm is
 * sensitive to roundoff, you might need to locally disable '-fp-model fast' by inserting this bit
 * of code at the start of each function that is sensitive:
 *
 * \code
 * int foo(int b) {
 *
 * #if defined __INTEL_COMPILER
 * #  pragma float_control(precise, on, push)
 * #endif
 *
 *   // code goes here
 *
 *   return value;
 * }
 * \endcode
 *
 * \sa [What Every Computer Scientist Should Know About Floating-Point
 * Arithmetic](https://docs.oracle.com/cd/E19957-01/806-3568/ncg_goldberg.html)
 */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_dsxx_FMA_hh
#define rtt_dsxx_FMA_hh

#include "ds++/config.h"
#include <cmath>

/*
 * Override - always use software fma or alwasy use hardware fma
 */
//#define FMA_ALWAYS_SOFTWARE 1
//#define FMA_NEVER_SOFTWARE 1

#if defined(FMA_ALWAYS_SOFTWARE) && defined(FMA_NEVER_SOFTWARE)
#error "Both FMA_ALWAYS_SOFTWARE and FMA_NEVER_SOFTWARE should not be set."
#endif

/*!
 * Normally \c FMA_FIND_DIFFS is \b NOT defined.  If defined then all FMA calls will be processed by
 * fma_with_diagnostics and will throw an exception if roundoff between the two computations is
 * significant.
 */
// #define FMA_FIND_DIFFS 1

#ifdef FMA_FIND_DIFFS

#include "ds++/Soft_Equivalence.hh"
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>

namespace rtt_dsxx {

//------------------------------------------------------------------------------------------------//
/*!
 * \brief fma_with_diagnostics
 *
 * \param[in] a Used when computing a*b + c
 * \param[in] b Used when computing a*b + c
 * \param[in] c Used when computing a*b + c
 * \param[in] file Filename where calculation occurs
 * \param[in] line Line number where calculation occurs
 * \param[in] abort_on_fail  Is failure fatal?
 * \return (a*b + c) if roundoff is small, otherwise throw an exception
 *
 * If \c abort_on_fail is \c true:
 * - Throw a C++ exception if the two values are sufficiently different
 * - Otherwise return the 'fast' result.
 *
 * If \c abort_on_fail is \c false:
 * - Print a message with details (don't throw) if the two values are
 *   sufficiently different.
 * - Return the 'accurate' result.
 */
template <typename T> // float, double, long double
inline T fma_with_diagnostics(T const a, T const b, T const c, std::string const &file,
                              uint32_t const line, bool const abort_on_fail) {
  T const accurate = fma(a, b, c);
  T const fast = ((a) * (b) + c);
  // For T=double, tol = 1.e+7 * 1.e-19 = 1.e-12
  T const tol =
      std::pow(10, (std::numeric_limits<T>::max_digits10 / 3)) * std::numeric_limits<T>::epsilon();
  std::ostringstream msg;
  msg << "FMA accuracy error in " << file << " (" << line << ") "
      << std::setprecision(std::numeric_limits<T>::max_digits10) << "\n  accurate = " << accurate
      << "\n  fast = " << fast << "\n  diff = " << fabs(accurate - fast);
  if (std::fabs(accurate + fast) < std::numeric_limits<T>::epsilon()) {
    // solution is zero, rdiff is meaningless
    msg << "\n  rdiff = NA";
  } else {
    msg << "\n  rdiff = "
        << 2.0 * fabs(accurate - fast) / fabs(accurate + fast + std::numeric_limits<T>::epsilon());
  }
  msg << "\n  a   = " << a << "\n  b   = " << b << "\n  c   = " << c << std::endl;
  if (abort_on_fail) {
    Insist(rtt_dsxx::soft_equiv(accurate, fast, tol), msg.str().c_str());
  } else {
    if (!rtt_dsxx::soft_equiv(accurate, fast, tol)) {
      std::cout << msg.str();
    }
  }

  return abort_on_fail ? fast : accurate;
}

} // namespace rtt_dsxx

//------------------------------------------------------------------------------------------------//
/*!
 * If \c FMA_FIND_DIFFS is defined, then hijack the whole FMA operation.  Instead, perform the
 * operation both ways and compare the result.  If the result is significantly different, throw an
 * exception.  This is useful for finding operations that are sensitive to roundoff.
 */
#define FMA(a, b, c) rtt_dsxx::fma_with_diagnostics((a), (b), (c), __FILE__, __LINE__, true)
#define FMA_ACCURATE(a, b, c)                                                                      \
  rtt_dsxx::fma_with_diagnostics((a), (b), (c), __FILE__, __LINE__, false)
#else

//------------------------------------------------------------------------------------------------//
//! HAVE_HARDWARE_FMA is set by config/platform_checks.cmake.
#ifdef HAVE_HARDWARE_FMA

//------------------------------------------------------------------------------------------------//
/*!
 * This will be the default for new machines. Both versions should be equally accurate, choose
 * whatever is faster. \c FP_FAST_FMA is set by the compiler
 *
 * \ref http://en.cppreference.com/w/c/numeric/math/fma
 *
 * Working on LANL systems, it seems that \c FP_FAST_FMA is not set by the Intel compiler, so Intel
 * will always pick the 2nd option.
 *
 * The intent of the \c FMA_ACCURATE version is to pick the more accurate multiply-plus-add
 * operation. When hardware FMA is available, we just choose the normal 'fast' version.
 */
#ifdef FP_FAST_FMA
#define FMA(a, b, c) fma((a), (b), (c))
#define FMA_ACCURATE(a, b, c) fma((a), (b), (c))
#else
#define FMA(a, b, c) ((a) * (b) + c)
#define FMA_ACCURATE(a, b, c) ((a) * (b) + c)
#endif

#else /* HAVE_HARDWARE_FMA is false */

//------------------------------------------------------------------------------------------------//
/*!
 * For older hardware that does \b not support FMA natively, provide a macro that will call an
 * infinite precision fma function so that numerical reproducibility is enhanced This allows
 * solutions generated on older machines to match solutions generated on newer machines.
 *
 * Operations that are known to be very sensitive to round-off error should use \c FMA_ACCURATE even
 * though it is much slower (2-3x slower).
 */
#define FMA(a, b, c) ((a) * (b) + c)
#define FMA_ACCURATE(a, b, c) fma((a), (b), (c))

#endif /* HAVE_HARDWARE_FMA */

#endif /* FMA_FIND_DIFFS */

//------------------------------------------------------------------------------------------------//
// Override
//------------------------------------------------------------------------------------------------//

#ifdef FMA_ALWAYS_SOFTWARE
#undef FMA
#define FMA(a, b, c) fma((a), (b), (c))
#undef FMA_ACCURATE
#define FMA_ACCURATE(a, b, c) fma((a), (b), (c))
#endif

#ifdef FMA_NEVER_SOFTWARE
#undef FMA
#define FMA(a, b, c) ((a) * (b) + c)
#undef FMA_ACCURATE
#define FMA_ACCURATE(a, b, c) ((a) * (b) + c)
#endif

#endif // rtt_dsxx_FMA_hh

//------------------------------------------------------------------------------------------------//
// end of ds++/FMA.hh
//------------------------------------------------------------------------------------------------//
