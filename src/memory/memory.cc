//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   memory/memory.cc
 * \author Kent G. Budge
 * \brief  memory diagnostic utilities
 * \note   Copyright (C) 2013-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "memory/memory.hh"
#include "ds++/Assert.hh"
#include "ds++/StackTrace.hh"
#include <iostream>
#include <limits>
#include <map>

#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wreserved-id-macro"
#pragma clang diagnostic ignored "-Wunused-macros"
#pragma clang diagnostic ignored "-Wunused-variable"
#endif

#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#endif

#ifndef _GLIBCXX_THROW
#define _GLIBCXX_THROW(except) throw(except)
#endif

namespace rtt_memory {
using namespace std;

static uint64_t total;
static uint64_t peak;
static uint64_t largest;

#if DRACO_DIAGNOSTICS & 2

// normally set in debugger to trigger a breakpoint
static uint64_t check_peak = numeric_limits<uint64_t>::max();

static uint64_t check_large = numeric_limits<uint64_t>::max();
// normally set in debugger to trigger a breakpoint

static uint64_t check_select_size = 504U;
// normally set in debugger to trigger a breakpoint

static uint64_t check_select_count = 0U;
// normally set in debugger to trigger a breakpoint

static unsigned dump_count = 0;
// indicates whether to dump the first few largest allocations and exit

static uint64_t report_threshold = numeric_limits<uint64_t>::max();

#endif

static bool is_active = false;

#if DRACO_DIAGNOSTICS & 2

struct alloc_t {
  std::size_t size; // size of allocation
  unsigned count;   // number of allocations of this size

  alloc_t() {}
  alloc_t(std::size_t my_size, unsigned my_count) : size(my_size), count(my_count) {}
};

struct Unsigned {
  unsigned value;

  Unsigned() : value(0) {}

  operator unsigned() const { return value; }
  unsigned &operator++() { return ++value; }
};

//================================================================================================//
/*!
 * \class memory_diagnostics
 * \brief
 *
 * We put the following in a wrapper so we can control destruction. We want to be sure is_active is
 * forced to be false once alloc_map is destroyed.
 */
//================================================================================================//
struct memory_diagnostics {
  map<void *, alloc_t> alloc_map;
  map<size_t, Unsigned> alloc_count;

  ~memory_diagnostics() { is_active = false; }
} st;

#endif // DRACO_DIAGNOSTICS & 2

//------------------------------------------------------------------------------------------------//
bool set_memory_checking(bool new_status) {
  bool Result = is_active;

#if DRACO_DIAGNOSTICS & 2
  total = 0;
  peak = 0;
  is_active = false;
  st.alloc_map.clear();
  st.alloc_count.clear();
#endif
  is_active = new_status;

  return Result;
}

//----------------------------------------------------------------------------------------//
#if DRACO_DIAGNOSTICS & 2
void set_report_threshold(uint64_t threshold) { report_threshold = threshold; }
#else
void set_report_threshold(uint64_t /*unused*/) {}
#endif

//------------------------------------------------------------------------------------------------//
uint64_t total_allocation() { return total; }

//------------------------------------------------------------------------------------------------//
uint64_t peak_allocation() { return peak; }

//------------------------------------------------------------------------------------------------//
uint64_t largest_allocation() { return largest; }

//------------------------------------------------------------------------------------------------//
/*! Print a report on possible leaks.
 *
 * This function prints a report in a human-friendly format on possible memory leaks.
 */
void report_leaks(ostream &out) {
  if (is_active) {
#if DRACO_DIAGNOSTICS & 2
    if (st.alloc_map.size() == 0) {
      out << "No indications of leaks" << endl;
    } else {
      map<void *, alloc_t>::const_iterator i;
      for (i = st.alloc_map.begin(); i != st.alloc_map.end(); ++i) {
        out << i->second.size << " bytes allocated at address " << i->first << " as allocation "
            << i->second.count << " of this size" << endl;
      }
    }
#else
    out << "No leak report available." << endl;
#endif
  }
}

//----------------------------------------------------------------------------------------//
// uint64_t set_check_peak(uint64_t new_peak) {
//   uint64_t Result = check_peak;
//   check_peak = new_peak;
//   return Result;
// }

//----------------------------------------------------------------------------------------//
// void set_dump_and_exit(unsigned new_dump_count) { dump_count = new_dump_count; }

} // end namespace rtt_memory

using namespace rtt_memory;

#if DRACO_DIAGNOSTICS & 2

//------------------------------------------------------------------------------------------------//
/*! Allocate memory with diagnostics.
 *
 * This version of operator new overrides the library default and allows us to track how memory is
 * being used while debugging. Since this introduces considerable overhead, it should not be used
 * for production builds.
 */
void *operator new(std::size_t n) _GLIBCXX_THROW(std::bad_alloc) {
  void *Result = malloc(n);

  // if malloc failed, then we need to deal with it.
  if (Result == 0) {
    // Store the global new handler
    // http://codereview.stackexchange.com/questions/7216/custom-operator-new-and-operator-delete
    bool failwithstacktrace(true);
    if (failwithstacktrace) {
      std::set_new_handler(rtt_memory::out_of_memory_handler);
      rtt_memory::out_of_memory_handler();
    }
    new_handler global_handler = set_new_handler(0);
    set_new_handler(global_handler);
    if (global_handler)
      global_handler();
    else
      throw bad_alloc();
  }

  // If malloc was successful, do the book keeping and return the pointer.
  if (is_active) {
    is_active = false;
    total += n;
    // Don't use max() here; doing it with if statement allows programmers to set a breakpoint here
    // to find high water marks of memory usage.
    if (total > peak) {
      peak = total;
      if (peak >= check_peak) {
        // This is where a programmer should set his breakpoint if he wishes to pause execution when
        // total memory exceeds the check_peak value (which the programmer typically also sets in
        // the debugger).
        cout << "Reached check peak value" << endl;
        if (dump_count > 0) {
          map<size_t, size_t> alloc_sizes;
          for (auto v : st.alloc_map) {
            if (alloc_sizes.find(v.second.size) == alloc_sizes.end()) {
              alloc_sizes[v.second.size] = 1;
            } else {
              alloc_sizes[v.second.size]++;
            }
          }
          auto iter = alloc_sizes.rbegin();
          if (alloc_sizes.find(n) == alloc_sizes.end()) {
            alloc_sizes[n] = 1;
          } else {
            alloc_sizes[n]++;
          }

          for (unsigned i = 0; i < dump_count; ++i) {
            cout << (*iter).first << ' ' << (*iter).second << endl;
            ++iter;
          }
          exit(EXIT_SUCCESS);
        }
      }
    }
    if (n >= check_large) {
      // This is where a programmer should set his breakpoint if he wishes to pause execution when a
      // memory allocation is requested that is larger than the check_large value (which the
      // programmer typically also sets in the debugger).
      cout << "Allocated check large value" << endl;
    }
    if (n > largest) {
      // Track the size of the largest single memory allocation.
      largest = n;
    }
    unsigned count = ++st.alloc_count[n];
    st.alloc_map[Result] = alloc_t(n, count);
    if (n == check_select_size && count == check_select_count) {
      // This is where the programmer should set his breakpoint if he wishes to pause execution on
      // the check_select_count'th instance of requesting an allocation of size check_select_size
      // (which the programmer typically also set in the debugger.) This is typically done to narrow
      // in on a potential memory leak, by identifying exactly which allocation is being leaked by
      // looking at the allocation map (st.alloc_map) to see the size and instance.
      cout << "Reached check select allocation" << endl;
    }
    if (n > report_threshold) {
      cout << "Memory allocation of size " << n << " was made. Total now " << total << '.' << endl;
    }
    is_active = true;
  }
  return Result;
}

//------------------------------------------------------------------------------------------------//
/*! Deallocate memory with diagnostics
 *
 * This is the operator delete override to go with the operator new override above.
 */
void operator delete(void *ptr) throw() {
  free(ptr);
  if (is_active) {
    map<void *, alloc_t>::iterator i = st.alloc_map.find(ptr);
    if (i != st.alloc_map.end()) {
      size_t const n = i->second.size;
      total -= n;
      if (n >= check_large) {
        // This is where the programmer should set his breakpoint if he wishes to pause execution
        // when an allocation larger than check_large is deallocated. check_large is typically also
        // set in the debugger by the programmer.
        cout << "Deallocated check large value" << endl;
      }
      if (n > report_threshold) {
        cout << "Memory allocation of size " << n << " was freed. Total now " << total << '.'
             << endl;
      }
      is_active = false;
      st.alloc_map.erase(i);
      is_active = true;
    }
  }
}

//------------------------------------------------------------------------------------------------//
/*! Deallocate memory with diagnostics
 *
 * C++14 introduces operator delete with a size_t argument, used in place of the unsized operator
 * delete when the size of the allocation can be deduced as a hint to the memory manager. For now,
 * we ignore the hint.
 *
 * Since C++14 does not mandate when the compiler calls this version, it is not possible to
 * guarantee coverage on all platforms.
 */
void operator delete(void *ptr, size_t) throw() { operator delete(ptr); }
#endif

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Provide a special action when an out-of-memory condition is encountered.
 *
 * The usual notion is that if new operator cannot allocate dynamic memory of the requested size,
 * then it should throw an exception of type std::bad_alloc.
 *
 * If std::bad_alloc is about to be thrown because new is unable to allocate enough memory, a
 * user-defined function can be called to provide diagnostic information.  This function must be
 * registered in the program.
 *
 * Example:
 *
 * \code
 * #include <cstdlib>
 * int main()
 * {
 *    // set the new handler.
 * #if DRACO_DIAGNOSTICS & 2
 *    std::set_new_handler(out_of_memory_handler);
 * #endif
 *    // invalid memory request
 *    int *pBigArray = new int[1000000000000L];
 *    return 0;
 * }
 * \endcode
 *
 * \bug untested
 */
void rtt_memory::out_of_memory_handler() {
  std::set_new_handler(nullptr);
  std::cerr << "Unable to allocate requested memory.\n" << rtt_dsxx::print_stacktrace("bad_alloc");
}

#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic pop
#endif

#ifdef __clang__
#pragma clang diagnostic pop
#endif

//------------------------------------------------------------------------------------------------//
// end of memory.cc
//------------------------------------------------------------------------------------------------//
