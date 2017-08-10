//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ds++/StackTrace.cc
 * \author Kelly Thompson
 * \date   Friday, Dec 20, 2013, 10:15 am
 * \brief  Linux/X86 implementation of stack trace functions.
 * \note   Copyright (C) 2016-2017 Los Alamos National Security, LLC.
 *         All rights reserved. */
//---------------------------------------------------------------------------//

#include "StackTrace.hh"
#include <iostream>
#include <sstream>

//---------------------------------------------------------------------------//
// Stack trace feature is only available on Unix-based systems when
// compiled with Intel or GNU C++.
//---------------------------------------------------------------------------//
#ifdef UNIX

#ifndef draco_isPGI
#include <cxxabi.h> // abi::__cxa_demangle
#endif
#include <execinfo.h> // backtrace
#include <stdio.h>    // snprintf
#include <stdlib.h>
#include <string.h>
#include <ucontext.h>
#include <unistd.h> // readlink

//---------------------------------------------------------------------------//
// Helper functions
std::string st_to_string(int const i) {
  std::ostringstream msg;
  msg << i;
  return msg.str();
}

//---------------------------------------------------------------------------//
// Print a demangled stack backtrace of the caller function.
std::string rtt_dsxx::print_stacktrace(std::string const &error_name) {
  // store/build the message here.  At the end of the function we return
  // msg.str().
  std::ostringstream msg;

  // Get our PID and build the name of the link in /proc
  pid_t const pid = getpid();

  // Build linkname
  std::string const linkname =
      std::string("/proc/") + st_to_string(pid) + std::string("/exe");

  // Now read the symbolic link (process name)
  unsigned const buf_size(512);
  char buf[buf_size];
  auto ret = readlink(linkname.c_str(), buf, buf_size);
  buf[ret] = 0;
  std::string process_name(buf);

  // retrieve current stack addresses
  int const max_frames = 64;
  void *addrlist[max_frames];
  int stack_depth = backtrace(addrlist, sizeof(addrlist) / sizeof(void *));

  // Print a header for the stack trace
  msg << "\nStack trace:"
      << "\n  Signaling error: " << error_name
      << "\n  Process        : " << process_name
      << "\n  PID            : " << pid
      << "\n  Stack depth    : " << stack_depth << " (showing "
      << stack_depth - 3 << ")"
      << "\n\n";

  // If there is no stack information, we are done.
  if (stack_depth == 0) {
    return msg.str();
  }

  // resolve addresses into strings containing "filename(function+address)",
  // this array must be free()-ed
  char **symbollist = backtrace_symbols(addrlist, stack_depth);

  // allocate string which will be filled with the demangled function name
  size_t funcnamesize = 256;
  char *funcname = (char *)malloc(funcnamesize);

  // msg << "\nRAW format:" << std::endl;
  // for( int i=0; i<stack_depth; ++i )
  // {
  //     msg << "  " << symbollist[i] << std::endl;
  // }
  // msg << "\nDemangled format:" << std::endl;

  // iterate over the returned symbol lines. skip first two,
  // (addresses of this function and handler)
  for (int i = 1; i < stack_depth - 2; i++) {
    char *begin_name = 0, *begin_offset = 0, *end_offset = 0;

    // find parentheses and +address offset surrounding the mangled name:
    // ./module(function+0x15c) [0x8048a6d]
    for (char *p = symbollist[i]; *p; ++p) {
      if (*p == '(')
        begin_name = p;
      else if (*p == '+')
        begin_offset = p;
      else if (*p == ')' && begin_offset) {
        end_offset = p;
        break;
      }
    }

    if (begin_name && begin_offset && end_offset && begin_name < begin_offset) {
      *begin_name++ = '\0';
      *begin_offset++ = '\0';
      *end_offset = '\0';
      char *location = end_offset + 1;

      // mangled name is now in [begin_name, begin_offset) and caller
      // offset in [begin_offset, end_offset). now apply
      // __cxa_demangle():

      int status(1); // assume failure
      char *ret = NULL;
#ifndef draco_isPGI
      ret = abi::__cxa_demangle(begin_name, funcname, &funcnamesize, &status);
#endif
      if (status == 0) {
        funcname = ret; // use possibly realloc()-ed string
        msg << "  " << symbollist[i] << " : " << funcname << "()+"
            << begin_offset << location << "\n";
      } else {
        // demangling failed. Output function name as a C function with
        // no arguments.
        msg << "  " << symbollist[i] << " : " << begin_name << "()+"
            << begin_offset << "\n";
      }
    } else {
      // couldn't parse the line? print the whole line.
      msg << "  " << symbollist[i] << " : ??\n";
    }
  }

  free(funcname);
  free(symbollist);

#ifdef draco_isPGI
  msg << "\n==> Draco's StackTrace feature is not currently implemented for "
         "PGI."
      << "\n    The StackTrace is known to work under Intel or GCC compilers."
      << std::endl;
#else
  msg << "\n==> Try to run 'addr2line -e " << process_name << " 0x99999' "
      << "\n    to find where each part of the stack relates to your source "
         "code."
      << "\n    Replace the 0x99999 with the actual address from the stack "
         "trace above."
      << std::endl;
#endif
  return msg.str();
}

#endif // UNIX

//---------------------------------------------------------------------------//
// Stack trace feature is also available on Win32-based systems when
// compiled with Visual Studio.
//---------------------------------------------------------------------------//
#ifdef WIN32

//---------------------------------------------------------------------------//
// Print a demangled stack backtrace of the caller function.
std::string rtt_dsxx::print_stacktrace(std::string const &error_name) {
  // store/build the message here.  At the end of the function we return
  // msg.str().
  std::ostringstream msg;

  int pid(0);
  std::string process_name("unknown");
  int stack_depth(3);

  // Print a header for the stack trace
  msg << "\nStack trace:"
      << "\n  Signaling error: " << error_name
      << "\n  Process        : " << process_name
      << "\n  PID            : " << pid
      << "\n  Stack depth    : " << stack_depth << " (showing "
      << stack_depth - 3 << ")"
      << "\n\n";

  msg << "\n==> Draco's StackTrace feature is not currently implemented for "
         "Win32."
      << "\n    The StackTrace is known to work under Intel or GCC compilers "
         "on Linux."
      << std::endl;

  return msg.str();
}

#endif // WIN32

//---------------------------------------------------------------------------//
// end of ds++/StackTrace.cc
//---------------------------------------------------------------------------//
