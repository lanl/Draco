//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   ds++/bin/testCrashHandler.cc
 * \author Kelly Thomposn
 * \date   Wed Nov 06 2013
 * \brief  Windows/x86 crash handler tests
 *
 * Copyright (C) 2014-2022 Triad National Security, LLC., All rights reserved.
 *
 * See ds++/COPYING file for more copyright information.  This code is based
 * substantially on http://www.codeproject.com/Articles/207464/Exception-Handling-in-Visual-Cplusplus
 */
//------------------------------------------------------------------------------------------------//

#include "ds++/Assert.hh"
#include "ds++/Release.hh"
#include "ds++/fpe_trap.hh"
#include <float.h> // _clearfp
#include <iostream>
#include <signal.h> // SIGABRT

#if defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-value"
#pragma clang diagnostic ignored "-Winfinite-recursion"
#pragma clang diagnostic ignored "-Wformat-security"
#endif

//------------------------------------------------------------------------------------------------//
void sigfpe_test() {
  // Code taken from http://www.devx.com/cplus/Article/34993/1954

  //Set the x86 floating-point control word according to what exceptions you want to trap.
  _clearfp(); //Always call _clearfp before setting the control word
  // Because the second parameter in the following call is 0, it only returns the floating-point
  // control word
  unsigned int cw;
  _controlfp_s(&cw, 0, 0); //Get the default control word
  // Set the exception masks off for exceptions that you want to trap.  When a mask bit is set, the
  // corresponding floating-point exception is blocked from being generating.
  cw &= ~(EM_OVERFLOW | EM_UNDERFLOW | EM_ZERODIVIDE | EM_DENORMAL | EM_INVALID);
  // For any bit in the second parameter (mask) that is 1, the corresponding bit in the first
  // parameter is used to update the control word.
  unsigned int cwOriginal;
  _controlfp_s(&cwOriginal, cw, MCW_EM); //Set it.

  // Divide by zero
  float a = 1;
  float b = 0;
  float c = a / b;
  c;
}

#define BIG_NUMBER 0x1fffffff
#pragma warning(disable : 4717) // avoid C4717 warning
int RecurseAlloc() {
  int *pi = new int[BIG_NUMBER];
  pi;
  RecurseAlloc();
  return 0;
}

class CDerived;
class CBase {
public:
  CBase(CDerived *derived) : m_pDerived(derived){};
  ~CBase();
  virtual void function(void) = 0;

  CDerived *m_pDerived;
};

#pragma warning(disable : 4355)
class CDerived : public CBase {
public:
  CDerived() : CBase(this){}; // C4355
  virtual void function(void){};
};

CBase::~CBase() { m_pDerived->function(); }

int main(int argc, char *argv[]) {
  // print a banner
  std::cout << rtt_dsxx::release() << std::endl;

  rtt_dsxx::CCrashHandler ch;
  ch.SetProcessExceptionHandlers();
  ch.SetThreadExceptionHandlers();

  int ExceptionType = 0;
  if (argc == 2) {
    ExceptionType = atoi(argv[1]);
  } else {
    std::cout << "Choose an exception type:\n"
              << "0 - SEH exception\n"
              << "1 - terminate\n"
              << "2 - unexpected\n"
              << "3 - pure virtual method call\n"
              << "4 - invalid parameter\n"
              << "5 - new operator fault\n"
              << "6 - SIGABRT\n"
              << "7 - SIGFPE\n"
              << "8 - SIGILL\n"
              << "9 - SIGINT\n"
              << "10 - SIGSEGV\n"
              << "11 - SIGTERM\n"
              << "12 - RaiseException\n"
              << "13 - throw C++ typed exception\n"
              << "Your choice >  ";
    std::cin >> ExceptionType; // scanf_s("%d", &ExceptionType);
  }

  std::cout << "Attempting to force ExceptionType = " << ExceptionType;

  switch (ExceptionType) {
  case 0: // SEH
  {
    // Access violation
    std::cout << " (SEH: Access violation)" << std::endl;
    int *p = 0;
#pragma warning(disable : 6011) // warning C6011: Dereferencing NULL pointer 'p'
    *p = 0;
#pragma warning(default : 6011)
    break;
  }
  case 1: {
    // Call terminate
    std::cout << " (terminate)" << std::endl;
    terminate();
    break;
  }
  case 2: {
    // Call unexpected
    std::cout << " (unexpected)" << std::endl;
    unexpected();
    break;
  }
  case 3: {
    // pure virtual method call
    std::cout << " (pure virtual method call)" << std::endl;
    CDerived derived;
    break;
  }
  case 4: {
    // invalid parameter
    std::cout << " (invalid parameter)" << std::endl;
    char *formatString;
    // Call printf_s with invalid parameters.
    formatString = NULL;
#pragma warning(                                                                                   \
    disable : 6387) // warning C6387: 'argument 1' might be '0': this does not adhere to the       \
                    // specification for the function 'printf'
    printf(formatString);
#pragma warning(default : 6387)
    break;
  }
  case 5: {
    // new operator fault: Cause memory allocation error
    std::cout << " (new operator fault)" << std::endl;
    RecurseAlloc();
    break;
  }
  case 6: {
    // Call abort (SIGABRT)
    std::cout << " (SIGABRT)" << std::endl;
    abort();
    break;
  }
  case 7: {
    // floating point exception ( /fp:except compiler option)
    std::cout << " (SIGFPE)" << std::endl;
    sigfpe_test();
    break;
  }
  case 8: {
    std::cout << " (SIGILL)" << std::endl;
    raise(SIGILL);
    break;
  }
  case 9: {
    std::cout << " (SIGINT)" << std::endl;
    raise(SIGINT);
    break;
  }
  case 10: {
    std::cout << " (SIGSEGV)" << std::endl;
    raise(SIGSEGV);
    break;
  }
  case 11: // SIGTERM
  {
    std::cout << " (SIGTERM)" << std::endl;
    raise(SIGTERM);
    break;
  }
  case 12: {
    // Raise noncontinuable software exception
    std::cout << " (Raise noncontinuable software exception)" << std::endl;
    RaiseException(123, EXCEPTION_NONCONTINUABLE, 0, NULL);
    break;
  }
  case 13: // throw
  {
    // Throw typed C++ exception.
    std::cout << " (Throw typed C++ exception)" << std::endl;
    throw 13;
    break;
  }
  default: {
    std::cout << " (Unknown exception type specified)" << std::endl;
    break;
  }
  }
  return 0;
}

#if defined(__clang__)
#pragma clang diagnostic pop
#endif

//------------------------------------------------------------------------------------------------//
// end of ds++/bin/testCrashHandler.cc
//------------------------------------------------------------------------------------------------//
