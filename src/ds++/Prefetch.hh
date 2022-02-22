//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   ds++/SystemCall.hh
 * \brief  Wrapper for explicit prefetch commands. Hides system dependence.
 * \note   Copyright (C) 2018-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_dsxx_Prefetch_hh
#define rtt_dsxx_Prefetch_hh

#include "ds++/config.h"

namespace rtt_dsxx {

//================================================================================================//
// Provides a common interface across architectures for explicit prefetches.
//================================================================================================//

/*! Prefetch a cache line.
 *
 * A call to this macro should be regarded as a hint to the processor that a particular cache line
 * is going to be read or written in the near future.  If hardware supports prefetching, the cache
 * line will begin to be read into cache while the processor carries out other computations, thus
 * ensuring that memory bus cycles are not wasted.
 *
 * An example would be where a significant computation is being performed on each element of an
 * array in turn. There is potential speedup if a prefetch is executed on the next element of the
 * array while the current element of the array is being computed upon.
 *
 * Prefetching by hand is tricky to get right. The hardware routinely prefetches instructions but
 * will also sometimes prefetch data when the access pattern is predictable. You cannot easily alter
 * hardware prefetching. The compiler can attempt to determine where prefetching may be profitable
 * and insert prefetch instructions. This actually slows the code in enough cases (the compiler
 * guesses wrong) that Intel disables compiler prefetching by default. You should first try enabling
 * compiler prefetching and see if it help before putting a lot of effort into hand prefetching.
 *
 * All that said, we have evidence that inserting prefetches by hand can noticeably speed up the
 * code where the compiler option fails, with enough effort and experimentation.
 *
 * Things to keep in mind:
 *
 * 1. Only so many memory fetches can be queued up at one time -- around 24 on the Intel KNL and
 * Haswell architectures, for example. If you insert a prefetch when the fetch queue is already
 * full, the processor stalls until the current memory cycle finishes and opens a slot. This is
 * functionally harmless but does not help performance. You will want to analyze your code to
 * estimate how many prefetches to issue so as not to saturate the prefetch slots.
 *
 * 2. You can execute a prefetch on any address, whether you have access to that address or
 * not. Thus prefetching past the end of an array is functionally harmless: There is no risk of a
 * segmentation violation even if reading or writing that address would yield such a violation. This
 * simplifies coding prefetches in loops. However, prefetching past the end of the array may waste a
 * memory cycle and cache line.
 *
 * 3. There's no point prefetching a cache line you will be using immediately.  You win only if
 * there are other computations for the processor to carry out using data already in cache while the
 * cache line is prefetched. In other words, prefetching is a way to parallelize the computation
 * units of the processor with the memory bus; if there is no parallelism there to be exploited, you
 * obviously can't exploit the parallelism.
 *
 * 4. There's only so much room in the nearer caches. If your prefetch into L1 flushes a line you
 * needed for immediate computation, you lose. On the other hand, if you can pull a line into L2
 * while the processor is computing on data already in L1, you win. Try low temporality arguments
 * first.
 *
 * This is implemented as a macro rather than an inline function because the second and third
 * arguments to __builtin_prefetch must be constexpr unsigned integers and there is no way to
 * express this with an inline function.
 *
 * \param addr Address of the cache line you wish to prefetch.
 *
 * \param for_write If zero, data in the cache line is going to be read (or read-modify-written.)
 * This is the usual case and so this is the default.  If the value is 1, the cache line is going to
 * be written only, so there is no need to read it first.
 *
 * \param temporality If 0, the data is going to be read or written only once in the near future and
 * so need not be cached. (So what's the point? Beats me.) If 3, the data is going to be used many
 * times in the near future and so should be pulled into all levels of cache. 2 means limited use in
 * the near future, so pull only into L2 or below; 1 means quite limited use in the near future, so
 * pull only into L3 (if it exists) or L2 (if there is no L3). We have some experience suggesting 1
 * is a good value to try first, but the default on g++ for the underlying __builtin_prefetch is 3,
 * which we hesitate to override.
 */

#ifdef __GNUC__

#define prefetch_cache_line(addr, for_write, temporality)                                          \
  __builtin_prefetch(addr, for_write, temporality)

#else
// not __GNUC__

/* No prefetching available. */
#define prefetch_cache_line(addr, for_write, temporality)

#endif

//! Number of bytes (char) on a cache line.

#ifdef draco_isKNL
unsigned const CACHE_LINE_CHAR = 64U;
#else
unsigned const CACHE_LINE_CHAR = 32U; // correct for Haswell; assumed for other
#endif

//! Number of ints in a cache line.

unsigned const CACHE_LINE_INT = CACHE_LINE_CHAR / sizeof(int);

//! Number of long ints in a cache line.

unsigned const CACHE_LINE_LONG = CACHE_LINE_CHAR / sizeof(long);

//! Number of double precision values in a cache line.

unsigned const CACHE_LINE_DOUBLE = CACHE_LINE_CHAR / sizeof(double);

} // namespace rtt_dsxx

#endif // rtt_dsxx_Prefetch_hh

//------------------------------------------------------------------------------------------------//
// end of Prefetch.hh
//------------------------------------------------------------------------------------------------//
