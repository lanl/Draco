//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   memory/memory.hh
 * \author Kent G. Budge
 * \brief  Memory utilities for diagnostic purposes.
 * \note   Copyright (C) 2013-2022 Triad National Security, LLC., All rights reserved.
 *
 * The memory utilities were written to address a need to identify the memory "high-water mark" in a
 * call sequence. This was not available with the existing memory checking tools. Other capabilities
 * gradually accreted themselves to this set of utilities, such as leak characterization.
 */
//------------------------------------------------------------------------------------------------//

#ifndef memory_memory_hh
#define memory_memory_hh

#include "ds++/config.h" // defines DRACO_DIAGNOSTICS
#include <cstdint>       // cstdint not available on PGI
#include <cstdlib>
#include <iostream>
#include <limits>

namespace rtt_memory {

bool set_memory_checking(bool new_status);
uint64_t total_allocation();
uint64_t peak_allocation();
uint64_t largest_allocation();
void report_leaks(std::ostream &out);

//! Register rtt_dsxx::print_stacktrace() as the respose to std::bad_alloc.
void out_of_memory_handler();

// To set the peak at which to break
// uint64_t set_check_peak(uint64_t new_peak);
void set_dump_and_exit(unsigned /*unused*/);

//! To get a report on the console of all allocations over a threshold size.
void set_report_threshold(uint64_t /*unused*/ = std::numeric_limits<uint64_t>::max());

} // namespace rtt_memory

#endif

//------------------------------------------------------------------------------------------------//
// end of memory/memory.hh
//------------------------------------------------------------------------------------------------//
