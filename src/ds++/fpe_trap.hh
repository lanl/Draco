//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   ds++/fpe_trap.hh
 * \author Rob Lowrie, Kelly Thompson
 * \date   Thu Oct 13 16:36:09 2005
 * \brief  Contains functions in the fpe_trap namespace.
 * \note   Copyright (C) 2016-2020 Triad National Security, LLC. All rights reserved.
 *
 * Copyright (C) 1994-2001  K. Scott Hunziker.
 * Copyright (C) 1990-1994  The Boeing Company.
 *
 * See ds++/COPYING file for more copyright information.  This code is based substantially on
 * fpe/i686-pc-linux-gnu.c from algae-4.3.6, which is available at http://algae.sourceforge.net/.
 */
//------------------------------------------------------------------------------------------------//

#ifndef fpe_trap_hh
#define fpe_trap_hh

#include "ds++/config.h"
#include <iostream>
#include <string>

namespace rtt_dsxx {

//------------------------------------------------------------------------------------------------//
/*!
 * \class fpe_trap
 *
 * \brief Enable trapping of floating-point exceptions.
 *
 * The floating-point exception behavior is platform dependent. Nevertheless, the goal of this class
 * is to turn on trapping for the following exceptions:
 *
 * - Division by zero.
 * - Invalid operation; for example, the square-root of a negative number.
 * - Overflow.
 *
 * If a floating-point exception is detected, the code will abort using a mode triggered by the
 * value of abortWithInsist.
 * - If true, ds++'s Insist is called; that is, a C++ exception is thrown.
 * - If false, the default mechanism defined by the compiler will be used. For most modern compilers
 *   this results in a stack trace.
 *
 * Typically, an application calls this function once, before any floating-point operations are
 * performed (e.g.: \c wedgehog/Function_Interfaces.cc). Note that all program functionality then
 * traps floating-point exceptions, including in libraries. Currently, there is no way to turn
 * trapping off once it has been turned on.
 *
 * Note: By Draco coding convention, fpe_traps are enabled when
 * \code
 * DRACO_DIAGNOSTIC && 4 == true.
 * \endcode
 *
 * Useful links:
 * - http://stackoverflow.com/questions/77005/how-to-generate-a-stacktrace-when-my-gcc-c-app-crashes
 * - With C++11, some of these features can easily be implemented in code. See
 *   http://en.cppreference.com/w/cpp/numeric/fenv/FE_exceptions .
 */
class fpe_trap {
public:
  /*!
   * \brief Constructor
   * \param[in] abortWithInsist_in toggle the abort mode default is true to use ds++'s Insist macro.
   */
  explicit fpe_trap(bool const abortWithInsist_in = true)
      : abortWithInsist(abortWithInsist_in) { /* emtpy */
  }

  //! Enable trapping of fpe signals.
  bool enable();
  //! Disable trapping of fpe signals.
  void disable();
  //! Query if trapping of fpe signals is active.
  bool active() const { return fpeTrappingActive; }

private:
  bool fpeTrappingActive{false};
  bool abortWithInsist;
};
} // end namespace rtt_dsxx

//------------------------------------------------------------------------------------------------//
/* WINDOWS X86
 *
 * Provide additional crash handlers for Win32 platforms. This code provides the ability to catch
 * and handle the following exceptions:
 * - Structured Exception Handling (SEH): An exception code is provided and typically triggers
 *   Dr.Watson to start.
 * - Vectored Exception Handling (VEH): These exceptions are supported in WinXP and later.  VEH
 *   allows the system to watch or handle all SEH chained exceptions for an application.
 * - CRT Error Handling: The C run time libraries provide their own error handling mechanism to
 *   catch C++ exceptions like Terminate, Pure Call, New Operator Fault and Invalid Parameter.
 * - C++ Signal Handling: C++ provides a program interruption mechanism to catch SIGABRT, SIGFPE,
 *   SIGILL, SIGINT, SIGSEGV and SIGTERM.
 *
 * http://www.codeproject.com/Articles/207464/Exception-Handling-in-Visual-Cplusplus
 *
 * In the class declaration, you also can see several exception handler functions, such as
 * SehHandler(), TerminateHandler() and so on. Any of these exception handlers can be called when an
 * exception occurs. A handler function (optionally) retrieves exception information and invokes
 * crash minidump generation code, then it terminates process with TerminateProcess() function call.
 *
 * The GetExceptionPointers() static method is used to retrieve exception information.
 */
//------------------------------------------------------------------------------------------------//
#ifdef FPETRAP_WINDOWS_X86
#include <Windows.h> // EXCEPTION_POINTERS

namespace rtt_dsxx {
class CCrashHandler {
public:
  // Constructor
  CCrashHandler(){/*empty*/};

  // Destructor
  virtual ~CCrashHandler(){/*empty*/};

  //! Sets exception handlers that work on per-process basis
  void SetProcessExceptionHandlers();

  //! Installs C++ exception handlers that function on per-thread basis
  void SetThreadExceptionHandlers();

  //! Collects current process state.
  static void GetExceptionPointers(DWORD dwExceptionCode, EXCEPTION_POINTERS **pExceptionPointers);

  //! Action to perform when an exception is found.
  static void ActionOnException(std::string const &message, EXCEPTION_POINTERS *pExcPtrs);

  /*!
   * \brief This method creates minidump of the process
   * \param pExcPtrs Pointer to the EXCEPTION_POINTERS structure containing
   *        exception information.
   *
   * The method calls the MiniDumpWriteDump() function from Microsoft Debug Help Library to
   * generate a minidump file.
   */
  static void CreateMiniDump(EXCEPTION_POINTERS *pExcPtrs);

  /* Exception handler functions. */

  static LONG WINAPI SehHandler(PEXCEPTION_POINTERS pExceptionPtrs);
  static void __cdecl TerminateHandler();
  static void __cdecl UnexpectedHandler();

  static void __cdecl PureCallHandler();

  static void __cdecl InvalidParameterHandler(const wchar_t *expression, const wchar_t *function,
                                              const wchar_t *file, unsigned int line,
                                              uintptr_t pReserved);

  static int __cdecl NewHandler(size_t);

  static void SigabrtHandler(int);
  static void SigfpeHandler(int /*code*/, int subcode);
  static void SigintHandler(int);
  static void SigillHandler(int);
  static void SigsegvHandler(int);
  static void SigtermHandler(int);
};

} // end namespace rtt_dsxx

#endif // end FPETRAP_WINDOWS_X86

#endif // fpe_trap_hh

//------------------------------------------------------------------------------------------------//
// end of ds++/fpe_trap.hh
//------------------------------------------------------------------------------------------------//
