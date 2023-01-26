/*---------------------------------------------*-C-*----------------------------------------------*/
/* file  : c4_omp.h
 * author: Kelly Thompson
 * date  : Wednesday, Nov 23, 2011, 08:42 am
 * brief : wrapper for omp.h
 * note  : Copyright (C) 2011-2023 Triad National Security, LLC., All rights reserved. */
/*------------------------------------------------------------------------------------------------*/

#ifndef rtt_c4_c4_omp_h
#define rtt_c4_c4_omp_h

#include "c4/config.h"

#ifdef OPENMP_FOUND
#include <omp.h>
#endif

/*------------------------------------------------------------------------------------------------*/
/* Wrap OpenMP API calls */
/*------------------------------------------------------------------------------------------------*/
namespace rtt_c4 {

#ifdef OPENMP_FOUND

/*! Set the number of OpenMP threads
 *\param n: desired number of threads
 */
inline void set_omp_num_threads(int n) { omp_set_num_threads(n); }

/*! Get the number of OpenMP threads in current section
 *\return: current number of threads
 */
inline int get_omp_num_threads() { return omp_get_num_threads(); }

/*! Get the maximum number of OpenMP threads
 *\return: maximum number of threads
 */
inline int get_omp_max_threads() { return omp_get_max_threads(); }

/*! Get the thread number for this thread
 *\return the thread number [0..N)
 */
inline int get_omp_thread_num() { return omp_get_thread_num(); }

/*! Set the schedule type and (optional) chunk size
 *
 * From https://www.openmp.org/spec-html/5.0/openmpsu121.html:
 *
 * typedef enum omp_sched_t {
 * omp_sched_static = 0x1,
 * omp_sched_dynamic = 0x2,
 * omp_sched_guided = 0x3,
 * omp_sched_auto = 0x4,
 * omp_sched_monotonic = 0x80000000u
 * }
 *
 * "For the schedule kinds static, dynamic, and guided the chunk_size is set to the value of the
 * second argument, or to the default chunk_size if the value of the second argument is less than 
 * 1."
 *
 * \param[in] kind Valid schedule type
 * \param[in] chunk_size (Optional) chunk size for scheduling
 */
inline void set_omp_schedule(const omp_sched_t kind, const int chunk_size = 0) {
  omp_set_schedule(kind, chunk_size);
}

/*! Get the current schedule type and chunk size
 *
 * \param[out] kind Valid schedule type
 * \param[out] chunk_size (Optional) chunk size for scheduling
 */
inline void get_omp_schedule(omp_sched_t &kind, int &chunk_size) {
  omp_get_schedule(&kind, &chunk_size);
}

#else

/*! Bypass setting the number of OpenMP threads */
inline void set_omp_num_threads(int /*n*/) { return; }

/*! Bypass getting the current number of OpenMP threads
 *\return: 1 (always)
 */
inline int get_omp_num_threads() { return 1; }

/*! Bypass getting the maximum number of OpenMP threads
 *\return: 1 (always)
 */
inline int get_omp_max_threads() { return 1; }

/*! Bypass getting the thread number for this thread
 *\return 0 (always)
 */
inline int get_omp_thread_num() { return 0; }

/*! Bypass setting schedule for OpenMP loops
 */
inline void set_omp_schedule(const int /*kind*/, const int /*chunk_size*/ = 0) { return; }

/*! Bypass getting schedule for OpenMP loops
 */
inline void get_omp_schedule(int & /*kind*/, int & /*chunk_size*/) { return; }
#endif

} // namespace rtt_c4

#endif /* rtt_c4_c4_omp_h */

/*------------------------------------------------------------------------------------------------*/
/* end of c4_omp.h */
/*------------------------------------------------------------------------------------------------*/
