//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   ds++/StackTrace.hh
 * \author Kelly Thompson
 * \date   Friday, Dec 20, 2013, 09:47 am
 * \brief  Contains function to generate a stack trace on Linux
 * \note   Copyright (C) 2014-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef StackTrace_hh
#define StackTrace_hh

#include "ds++/config.h"
#include <string>

namespace rtt_dsxx {

//------------------------------------------------------------------------------------------------//
/*!
 * \brief print_stacktrace
 *
 * \param error_name A string that identifies why the stack trace is requested.
 * \return A multi-line message including the error_name and the stack trace.
 *
 * \note Consider replacing this with C++23's 'stacktrace' feature.
 *
 * A stack trace will look something like this:
 * \code
 Stack trace:
 Signaling error: myMessage
 Process        : /var/tmp/kgt/gcc-mpid-diag3/d/src/ds++/test/tstStackTrace
 PID            : 29732
 Stack depth    : 7 (showing 4)

 ./tstStackTrace : sr2(std::string&)()+0x41 [0x405671]
 ./tstStackTrace : sr1(std::string&)()+0x18 [0x405705]
 ./tstStackTrace : runtest(rtt_dsxx::UnitTest&)()+0x54 [0x40575c]
 ./tstStackTrace : main()+0x50
 * \code
 *
 * To debug a stack trace use objdump and c++filt (Intel/g++ only).
 * \code
 * objdump -dS tstStackTrace | c++filt > tstStackTrace.dmp
 * \endcode
 * View the dump file, searching for the the offset value (0x50).
 *
 * Or, you can use addr2line:
 * \code
 * addr2line -e ./tstStackTrace 0x405671
 * /home/kellyt/draco/src/ds++/test/tstStackTrace.cc:29
 * addr2line -e ./tstStackTrace 0x405705
 * /home/kellyt/draco/src/ds++/test/tstStackTrace.cc:36
 * \endcode
 *
 * Useful links:
 * - http://stackoverflow.com/questions/77005/how-to-generate-a-stacktrace-when-my-gcc-c-app-crashes
 * - http://stackwalker.codeplex.com
 */
std::string print_stacktrace(std::string const &error_name);

} // end namespace rtt_dsxx

#endif // StackTrace_hh

//------------------------------------------------------------------------------------------------//
// end of ds++/StackTrace.hh
//------------------------------------------------------------------------------------------------//
