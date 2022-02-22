//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   c4/C4_Serial.cc
 * \author Thomas M. Evans
 * \date   Mon Mar 25 17:06:25 2002
 * \brief  Implementation of C4 serial option.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "C4_Functions.hh"
#include "ds++/SystemCall.hh"
#include <chrono>
#include <cstdlib>
#include <ctime>

#ifdef C4_SCALAR

namespace rtt_c4 {

//------------------------------------------------------------------------------------------------//
// Any source rank
//------------------------------------------------------------------------------------------------//

const int any_source = -1;

//------------------------------------------------------------------------------------------------//
// Null source/destination rank
//------------------------------------------------------------------------------------------------//

const int proc_null = -2;

//------------------------------------------------------------------------------------------------//
// SETUP FUNCTIONS
//------------------------------------------------------------------------------------------------//

int initialize(int & /* argc */, char **& /*argv */, int /*required*/) { return 0; }

//------------------------------------------------------------------------------------------------//
void finalize() {}

//------------------------------------------------------------------------------------------------//
void type_free(C4_Datatype & /*old_type*/) {}

//------------------------------------------------------------------------------------------------//
void free_inherited_comm() {}

//------------------------------------------------------------------------------------------------//
// QUERY FUNCTIONS
//------------------------------------------------------------------------------------------------//

int node() { return 0; }
uint32_t rank() { return 0; }

//------------------------------------------------------------------------------------------------//
int nodes() { return 1; }
uint32_t nranks() { return 1; }

//------------------------------------------------------------------------------------------------//
// BARRIER FUNCTIONS
//------------------------------------------------------------------------------------------------//

void global_barrier() { /* empty */
}

//------------------------------------------------------------------------------------------------//
// TIMING FUNCTIONS
//------------------------------------------------------------------------------------------------//

#if defined(_MSC_VER)
double wall_clock_time(DRACO_TIME_TYPE &now) {
  using namespace std::chrono;
  now = high_resolution_clock::now();
  high_resolution_clock::duration t0 = now.time_since_epoch();
  return (static_cast<double>(t0.count()) * high_resolution_clock::period::num) /
         high_resolution_clock::period::den;
}
double wall_clock_time() {
  DRACO_TIME_TYPE now;
  return wall_clock_time(now);
}
#else
/* Linux */
double wall_clock_time() {
  DRACO_TIME_TYPE now;
  return static_cast<double>(times(&now)) / wall_clock_resolution();
}
double wall_clock_time(DRACO_TIME_TYPE &now) {
  return static_cast<double>(times(&now)) / wall_clock_resolution();
}
#endif

//------------------------------------------------------------------------------------------------//

double wall_clock_resolution() { return static_cast<double>(DRACO_CLOCKS_PER_SEC); }

//------------------------------------------------------------------------------------------------//
// PROBE/WAIT FUNCTIONS
//------------------------------------------------------------------------------------------------//

bool probe(int /* source */, int /* tag */, int & /* message_size */) { return false; }

void blocking_probe(int /* source */, int /* tag */, int & /* message_size */) {
  Insist(false, "no messages expected in serial programs!");
}

void wait_all(unsigned /*count*/, C4_Req *const /*requests*/) { return; }

std::vector<int> wait_all_with_source(unsigned /*count*/, C4_Req *const /*requests*/) {
  return std::vector<int>();
}

unsigned wait_any(unsigned /*count*/, C4_Req * /*requests*/) {
  Insist(false, "no messages expected in serial programs!");
  return 0;
}

//------------------------------------------------------------------------------------------------//
// ABORT
//------------------------------------------------------------------------------------------------//

#ifdef _MSC_FULL_VER
// - 4702 unreachable code
#pragma warning(push)
#pragma warning(disable : 4702)
#endif

int abort(int error) {
  // call system exit
  std::abort();
  return error; // unreachable
}

#ifdef _MSC_FULL_VER
#pragma warning(pop)
#endif

//------------------------------------------------------------------------------------------------//
// isScalar
//------------------------------------------------------------------------------------------------//
bool isScalar() { return true; }
bool isMpiInit() { return false; }

//------------------------------------------------------------------------------------------------//
// get_processor_name
//------------------------------------------------------------------------------------------------//
std::string get_processor_name() {
  return rtt_dsxx::draco_gethostname();
  ;
}

} // end namespace rtt_c4

#endif // C4_SCALAR

//------------------------------------------------------------------------------------------------//
// end of C4_Serial.cc
//------------------------------------------------------------------------------------------------//
