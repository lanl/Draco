//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   ds++/StackTrace.cc
 * \author Kelly Thompson
 * \date   Friday, Dec 20, 2013, 10:15 am
 * \brief  Linux/X86 implementation of stack trace functions.
 * \note   Copyright (C) 2014-2021 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "StackTrace.hh"
#include <array>
#include <iostream>
#include <sstream>

//------------------------------------------------------------------------------------------------//
// Stack trace feature is only available on Unix-based systems when compiled with Intel or GNU C++.
//------------------------------------------------------------------------------------------------//
#ifdef UNIX

#ifndef draco_isPGI
#include <cxxabi.h> // abi::__cxa_demangle
#endif
#include <cstdio> // snprintf
#include <cstdlib>
#include <cstring>
#include <execinfo.h> // backtrace
#include <ucontext.h>
#include <unistd.h> // readlink

//------------------------------------------------------------------------------------------------//
// Helper functions
inline std::string st_to_string(int const i) {
  std::ostringstream msg;
  msg << i;
  return msg.str();
}

//------------------------------------------------------------------------------------------------//
// Print a demangled stack backtrace of the caller function.
std::string rtt_dsxx::print_stacktrace(std::string const &error_message) {
  // store/build the message here.  At the end of the function we return msg.str().
  std::ostringstream msg;

  // Get our PID and build the name of the link in /proc
  pid_t const pid = getpid();

  // Now read the symbolic link (process name)
  unsigned constexpr buf_size(512);
  std::array<char, buf_size> buf;
#ifdef APPLE
  int ret = -1; // This scheme won't work on OSX: no /proc fs
#else
  // Build linkname
  std::string const linkname = std::string("/proc/") + st_to_string(pid) + std::string("/exe");

  auto ret = readlink(linkname.c_str(), buf.data(), buf_size);
#endif
  if (ret >= 0) /* readlink succeeded */
    buf[ret] = 0;
  std::string const process_name = ret < 0 ? "UNAVAILABLE" : std::string(buf.data());

  // retrieve current stack addresses
  int constexpr max_frames = 64;
  std::array<void *, max_frames> addrlist;
  uint32_t constexpr sizeofvoidptr = sizeof(void *);
  int const stack_depth = backtrace(addrlist.data(), sizeofvoidptr * addrlist.size());

  // Print a header for the stack trace
  msg << "\n"
      << error_message << "\nStack trace:"
      << "\n  Process        : " << process_name << "\n  PID            : " << pid
      << "\n  Stack depth    : " << stack_depth << " (showing " << stack_depth - 2 << ")"
      << "\n\n";

  // If there is no stack information, we are done.
  if (stack_depth == 0) {
    return msg.str();
  }

  // resolve addresses into strings containing "filename(function+address)", this array must be
  // free()-ed
  char **symbollist = backtrace_symbols(addrlist.data(), stack_depth);

  // allocate string which will be filled with the demangled function name
  size_t funcnamesize = 256;
  auto *funcname = static_cast<char *>(malloc(funcnamesize));

  // iterate over the returned symbol lines. skip first two, (addresses of this function and
  // handler)
  for (int i = 0; i < stack_depth - 2; i++) {
    char *begin_name = nullptr;
    char *begin_offset = nullptr;
    char *end_offset = nullptr;

    // find parentheses and +address offset surrounding the mangled name: ./module(function+0x15c)
    // [0x8048a6d]
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

      // mangled name is now in [begin_name, begin_offset) and caller offset in [begin_offset,
      // end_offset). now apply __cxa_demangle():

      int status(1); // assume failure
      char *ret01 = nullptr;
#ifndef draco_isPGI
      ret01 = abi::__cxa_demangle(begin_name, funcname, &funcnamesize, &status);
#endif
      if (status == 0) {
        funcname = ret01; // use possibly realloc()-ed string
        msg << "  " << symbollist[i] << " : " << funcname << "()+" << begin_offset << location
            << "\n";
      } else {
        // demangling failed. Output function name as a C function with no arguments.
        msg << "  " << symbollist[i] << " : " << begin_name << "()+" << begin_offset << "\n";
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
      << "\n    The StackTrace is known to work under Intel or GCC compilers." << std::endl;
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

//------------------------------------------------------------------------------------------------//
// Stack trace feature is also available on Win32-based systems when compiled with Visual Studio.
// ------------------------------------------------------------------------------------------------//
#ifdef _MSC_VER

//------------------------------------------------------------------------------------------------//
// Print a demangled stack backtrace of the caller function.
std::string rtt_dsxx::print_stacktrace(std::string const &error_message) {
  // store/build the message here.  At the end of the function we return msg.str().
  std::ostringstream msg;

  int pid(0);
  std::string process_name("unknown");
  int stack_depth(3);

  // Print a header for the stack trace
  msg << "\n"
      << error_message << "\nStack trace:"
      << "\n  Process        : " << process_name << "\n  PID            : " << pid
      << "\n  Stack depth    : " << stack_depth << " (showing " << stack_depth - 3 << ")"
      << "\n\n";

  msg << "\n==> Draco's StackTrace feature is not currently implemented for "
         "Win32."
      << "\n    The StackTrace is known to work under Intel or GCC compilers "
         "on Linux."
      << std::endl;

  return msg.str();
}

#endif // _MSC_VER

//------------------------------------------------------------------------------------------------//
// end of ds++/StackTrace.cc
//------------------------------------------------------------------------------------------------//
