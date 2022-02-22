//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   c4/Global_Timer.cc
 * \author Kent G. Budge
 * \date   Mon Mar 25 17:35:07 2002
 * \brief  Define methods of class Global_Timer, a POSIX standard timer.
 * \note   Copyright (C) 2013-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "Global_Timer.hh"

namespace rtt_c4 {
using namespace std;

bool Global_Timer::global_active_ = false;
map<string, Global_Timer::timer_entry> Global_Timer::active_list_;

//------------------------------------------------------------------------------------------------//
Global_Timer::Global_Timer(char const *name) : name_(name), active_(false) {
  Require(name != nullptr);
  timer_entry &entry = active_list_[name];
  active_ = entry.is_active;
  Check(entry.timer == nullptr); // Global_Timers must have unique names.
  entry.timer = this;

  Ensure(name == this->name());
}

//------------------------------------------------------------------------------------------------//
void Global_Timer::set_selected_activity(set<string> const &timer_list, bool const active) {
  if (rtt_c4::node() == 0) {
    cout << "***** Global timers selectively activated:" << endl;
    for (auto const &name : timer_list) {
      cout << " \"" << name << '\"';
      timer_entry &entry = active_list_[name];
      if (entry.timer == nullptr) {
        cout << " (DEFERRED)";
      }
      cout << endl;
    }
  }
  for (auto const &name : timer_list) {
    timer_entry &entry = active_list_[name];
    entry.is_active = active;
    if (entry.timer != nullptr) {
      entry.timer->set_activity(active);
    }
  }
}

//------------------------------------------------------------------------------------------------//
void Global_Timer::set_global_activity(bool const active) {
  global_active_ = active;
  std::string const status = active ? "ACTIVE" : "INACTIVE";
  if (rtt_c4::node() == 0)
    cout << "***** Global timers are now " << status << endl;
}

//------------------------------------------------------------------------------------------------//
void Global_Timer::reset_all() {
  if (rtt_c4::node() == 0)
    cout << "***** Resetting all global timers" << endl;
  for (auto const &i : active_list_) {
    timer_entry const entry = i.second;
    if ((entry.is_active || global_active_) && entry.timer) {
      entry.timer->reset();
    }
  }
}

//------------------------------------------------------------------------------------------------//
void Global_Timer::report_all(ostream &out) {
  if (rtt_c4::node() == 0) {

    std::string const divider(92U, '-');
    out << divider << "\nTiming report for all global timers:\n\n"
        << "            N                   user                   system                wall"
        << endl;
  }
  bool deferred_not_found = false;
  for (auto const &i : active_list_) {
    timer_entry const entry = i.second;
    if (entry.is_active || global_active_) {
      if (entry.timer != nullptr) {
        if (rtt_c4::node() == 0) {
          out << entry.timer->name() << endl;
        }
        entry.timer->printline_mean(out);
      } else {
        deferred_not_found = true;
      }
    }
  }
  if (rtt_c4::node() == 0) {
    if (!global_active_ && deferred_not_found) {
      out << "**** WARNING: DEFERRED timers not found:" << endl;
      for (auto const &i : active_list_) {
        timer_entry const entry = i.second;
        if (entry.is_active && entry.timer == nullptr) {
          out << i.first << endl;
        }
      }
      out << endl << "Perhaps you wanted one of these timers found but not active?" << endl;
      for (auto const &i : active_list_) {
        timer_entry const entry = i.second;
        if (!entry.is_active && entry.timer != nullptr) {
          out << i.first << endl;
        }
      }
    }
    out << string(92U, '-') << endl;
  }
}

} // end namespace rtt_c4

//------------------------------------------------------------------------------------------------//
// end of c4/Global_Timer.cc
//------------------------------------------------------------------------------------------------//
