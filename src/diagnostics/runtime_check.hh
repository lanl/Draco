//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   diagnostics/runtime_check.hh
 * \author Kent Grimmett Budge
 * \brief  Define runtime_check function
 * \note   Copyright (C) 2018-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef diagnostics_runtime_check_hh
#define diagnostics_runtime_check_hh

#include "ds++/config.h"

namespace rtt_diagnostics {

//------------------------------------------------------------------------------------------------//
/*! Parallel synchronous runtime check
*
* This function supports a synchronous check of a local condition on all processors. If any
* processor fails the check, all will throw a std::runtime_error with the prescribed message. This
* is much cleaner than throwing a std::runtime_error on only the processors that fail, which will
* not permit any kind of recovery and leaves it to MPI to notice that a process has died and kill
* all the others. (Or worse yet, hang waiting for a message that will never arrive.
*
* This function must be called synchronously on all processors. It is intended for handling of run
* time errors in production code such as i/o failures or numerical breakdown that are not
* particularly a result of a code defect.
*
* \param condition Condition to be tested locally. If this is true on all processors, the function
*                  returns normally. If this is false on any processor, all processors throw 
*                  \c std::runtime_error with a message based on:
*
* \param message Message describing the condition being tested. This will be incorporated in the
*                message text for the std::runtime_error thrown on all processors if any fail the 
*                condition being tested.
*/
void runtime_check(bool condition, char const *message) noexcept(false);

} // end namespace rtt_diagnostics

#endif // diagnostics_runtime_check_hh

//------------------------------------------------------------------------------------------------//
// end of diagnostics/runtime_check.hh
//------------------------------------------------------------------------------------------------//
