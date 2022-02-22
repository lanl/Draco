//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   diagnostics/Timing.cc
 * \author Thomas M. Evans
 * \date   Tue Dec 13 10:44:29 2005
 * \brief  Timing class member definitions.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "Timing.hh"
#include "ds++/Assert.hh"
#include "ds++/Soft_Equivalence.hh"

namespace rtt_diagnostics {

//------------------------------------------------------------------------------------------------//
// STATIC PUBLIC FUNCTIONAL INTERFACE
//------------------------------------------------------------------------------------------------//
/*!
 * \brief Add to a specified timer.
 *
 * This functions adds value to the timer with name key.  The first time this is called the value is
 * added to zero.  The timers are all static, so multiple calls to the same key will keep a running
 * tally.  To reset, call reset_timer().
 *
 * Calling this function adds the timer with name key to the map of timers.
 */
void Timing_Diagnostics::update_timer(const std::string &key, double value) {
  timers[key] += value;
}

//------------------------------------------------------------------------------------------------//
//! Return a vector of timer keys.
Timing_Diagnostics::Vec_Keys Timing_Diagnostics::timer_keys() {
  Vec_Keys keys(timers.size()); // keys
  auto v = keys.begin();        // iterators

  // add keys to the vector
  for (auto m = timers.begin(); m != timers.end(); m++, v++) {
    Check(v != keys.end());
    *v = m->first;
  }

  // return the vector
  return keys;
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Reset a timer to zero.
 *
 * Calling this function adds the timer with name key to the map of timers.
 */
void Timing_Diagnostics::reset_timer(const std::string &key) {
  timers[key] = 0.0;
  Ensure(rtt_dsxx::soft_equiv(timers[key], 0.0));
}

//------------------------------------------------------------------------------------------------//
//! Reset all timers in the map of timers to zero.
void Timing_Diagnostics::reset_timers() {
  // reset each timer
  for (auto &m : timers)
    m.second = 0.0;
}

//------------------------------------------------------------------------------------------------//
//! Removes a timer with name key from the map of timers.
void Timing_Diagnostics::delete_timer(const std::string &key) {
  timers.erase(key);
  Ensure(timers.count(key) == 0);
}

//------------------------------------------------------------------------------------------------//
//! Delete all timers from the map.
void Timing_Diagnostics::delete_timers() {
  // null map
  std::map<std::string, double> null;

  // swap it with timers
  timers.swap(null);
  Ensure(timers.empty());
}

//------------------------------------------------------------------------------------------------//
// PRIVATE STATIC CLASS-MEMBER DEFINITIONS
//------------------------------------------------------------------------------------------------//

std::map<std::string, double> Timing_Diagnostics::timers;

} // end namespace rtt_diagnostics

//------------------------------------------------------------------------------------------------//
// end of Timing.cc
//------------------------------------------------------------------------------------------------//
