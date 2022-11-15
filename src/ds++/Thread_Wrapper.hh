//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   ds++/Thread_Wrapper.hh
 * \author Tim Kelley
 * \date   Thursday, Oct 12, 2017, 10:50 am
 * \brief  Header file for Thread_Wrapper
 * \note   Copyright (C) 2017-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef THREAD_WRAPPER_HH
#define THREAD_WRAPPER_HH

#include "Assert.hh"
#include <thread>

namespace rtt_dsxx {

//================================================================================================//
/**\class Thread_Wrapper
 * \brief Joins underlying thread in dtor, avoiding untimely process termination.
 *
 * cf Scott Meyers, Effective Modern C++, Item #37
 */
//================================================================================================//
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
  Thread_Wrapper() : t_(std::thread()) {}

  Thread_Wrapper(Thread_Wrapper const &rhs) = delete;
  Thread_Wrapper(Thread_Wrapper const &&rhs) noexcept = delete;
  Thread_Wrapper &operator=(Thread_Wrapper const &rhs) = delete;
  Thread_Wrapper &operator=(Thread_Wrapper &&rhs) noexcept = delete;

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
  action const action_{action::join};
  std::thread t_;
};
} // namespace rtt_dsxx

#endif // include guard

//------------------------------------------------------------------------------------------------//
// end of ds++/Thread_Wrapper.hh
//------------------------------------------------------------------------------------------------//
