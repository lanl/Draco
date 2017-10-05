//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ds++/Thread_Wrapper.hh
 * \author Tim Kelley
 * \date   Mar 20, 2017
 * \brief  Header file for Thread_Wrapper
 * \note   Copyright (C) 2006-2017 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
//---------------------------------------------------------------------------//

#ifndef THREAD_WRAPPER_HH
#define THREAD_WRAPPER_HH

#include "Assert.hh"
#include <thread>

namespace rtt_dsxx {

//===========================================================================//
/**\class Thread_Wrapper
 * \brief Joins underlying thread in dtor, avoiding untimely process termination.
 *
 * cf Scott Meyers, Effective Modern C++, Item #37
 */
//===========================================================================//
class Thread_Wrapper {
public:
  /**\brief enum class to indicate action to take on termination. */
  enum class action {
    join,  /*!< roughly, wait for the thread to finish. */
    detach /*!< detach from the thread: i.e. leave it to continue running.*/
  };

  /**\brief accessor */
  std::thread &get() { return t_; }

  /**\brief Default ctor*/
  Thread_Wrapper() : action_(action::join), t_(std::thread()) {}

  /**\brief Move ctor (note: no copy ctor on std::threads. */
  Thread_Wrapper(std::thread &&t, action a) : action_(a), t_(std::move(t)) {}

  ~Thread_Wrapper() {
    if (t_.joinable()) {
      if (action::join == action_) {
        t_.join();
      } else {
        t_.detach();
      }
    }
    Check(!t_.joinable());
  } // dtor

private:
  action const action_;
  std::thread t_;
};
} // rtt_dsxx::

#endif // include guard

//---------------------------------------------------------------------------//
// end of ds++/Thread_Wrapper.hh
//---------------------------------------------------------------------------//
