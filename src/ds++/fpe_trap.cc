//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   ds++/fpe_trap.cc
 * \author Rob Lowrie, Kelly Thompson
 * \date   Thu Oct 13 16:52:05 2005
 * \brief  platform dependent implementation of fpe_trap functions.
 * \note   Copyright (C) 2013-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "fpe_trap.hh"
#include "Assert.hh"
#include "StackTrace.hh"
#include <sstream>

//------------------------------------------------------------------------------------------------//
// Linux_x86
//------------------------------------------------------------------------------------------------//
#ifdef FPETRAP_LINUX_X86

#include <cfenv>
#include <csignal>

/* Signal handler for floating point exceptions. */
extern "C" [[noreturn]] void catch_sigfpe(int sig, siginfo_t *psSiginfo, void * /*psContext*/) {
  // generate a message:
  std::string error_type;

  if (sig != SIGFPE) {
    error_type = "FATAL ERROR: Floating point exception problem.";
  } else {
    switch (psSiginfo->si_code) {
    case FPE_INTDIV:
      error_type = "FATAL ERROR (SIGNAL) Caught SIGFPE (Integer divide by zero)";
      break;
    case FPE_INTOVF:
      error_type = "FATAL ERROR (SIGNAL) Caught SIGFPE (Integer overflow)";
      break;
    case FPE_FLTDIV:
      error_type = "FATAL ERROR (SIGNAL) Caught SIGFPE (Floating point divide by zero)";
      break;
    case FPE_FLTOVF:
      error_type = "FATAL ERROR (SIGNAL) Caught SIGFPE (Floating point overflow)";
      break;
    case FPE_FLTUND:
      error_type = "FATAL ERROR (SIGNAL) Caught SIGFPE (Floating point underflow)";
      break;
    case FPE_FLTRES:
      error_type = "FATAL ERROR (SIGNAL) Caught SIGFPE (Floating point inexact result)";
      break;
    case FPE_FLTINV:
      error_type = "FATAL ERROR (SIGNAL) Caught SIGFPE (Invalid floating point "
                   "operation)";
      break;
    case FPE_FLTSUB:
      error_type = "FATAL ERROR (SIGNAL) Caught SIGFPE (Floating point "
                   "subscript out of range)";
      break;
    default:
      error_type = "FATAL ERROR (SIGNAL) Caught SIGFPE (Unknown floating point "
                   "exception)";
      break;
    }
  }
  Insist(false, rtt_dsxx::print_stacktrace(error_type));
}

//------------------------------------------------------------------------------------------------//
namespace rtt_dsxx {

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Enable trapping fpe signals.
 * \return \b true if trapping is enabled, \b false otherwise.
 *
 * A \b false return value is typically because the platform is not supported.
 */
bool fpe_trap::enable() {
  struct sigaction act {};

  // Choose to use Draco's DbC Insist.  If set to false, the compiler should print a stack trace
  // instead of the pretty print message defined above in catch_sigfpe.
  if (this->abortWithInsist)
    act.sa_sigaction = catch_sigfpe; // the signal handler

  sigemptyset(&(act.sa_mask)); // no other signals blocked
  act.sa_flags = SA_SIGINFO;   // want 3 args for handler

  // specify handler
  Insist(!sigaction(SIGFPE, &act, nullptr), "Unable to set floating point handler.");

  // The feenableexcept function is new for glibc 2.2.  See its description in the man page for
  // fenv(3).

  (void)feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW); // NOLINT [hicpp-signed-bitwise]

  // Toggle the state.
  fpeTrappingActive = true;

  // ... if you want to catch "everything" :
  //(void)feenableexcept(FE_ALL_EXCEPT);

  // Also available (taken from fenv.h and bits/fenv.h):
  // FE_INEXACT               inexact result
  // FE_DIVBYZERO             division by zero
  // FE_UNDERFLOW             result not representable due to underflow
  // FE_OVERFLOW              result not representable due to overflow
  // FE_INVALID               invalid operation

  // FE_ALL_EXCEPT            bitwise OR of all supported exceptions

  // The next macros are defined iff the appropriate rounding mode is supported by the
  // implementation.
  // FE_TONEAREST             round to nearest
  // FE_UPWARD                round toward +Inf
  // FE_DOWNWARD              round toward -Inf
  // FE_TOWARDZERO            round toward 0

  return fpeTrappingActive;
}

//------------------------------------------------------------------------------------------------//
//! Disable trapping of floating point errors.
void fpe_trap::disable() {
  (void)fedisableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW); // NOLINT [hicpp-signed-bitwise]
  fpeTrappingActive = false;
  return;
}

} // namespace rtt_dsxx

#endif // FPETRAP_LINUX_X86

//------------------------------------------------------------------------------------------------//
// OSF_ALPHA
//------------------------------------------------------------------------------------------------//
#ifdef FPETRAP_OSF_ALPHA

#include <machine/fpu.h>
#include <signal.h>

// Local functions

extern "C" {

/* Signal handler for floating point exceptions. */

static void catch_sigfpe(int sig) {
  std::string mesg = "Floating point exception";
  // decipher sig later...
  Insist(0, mesg);
}

} // end of namespace

namespace rtt_dsxx {

//------------------------------------------------------------------------------------------------//
//!  Enable trapping of floating point errors.
bool fpe_trap::enable() {
  unsigned long csr = ieee_get_fp_control();
  csr |= IEEE_TRAP_ENABLE_INV | IEEE_TRAP_ENABLE_DZE | IEEE_TRAP_ENABLE_OVF;
  ieee_set_fp_control(csr);

  if (abortWithInsist)
    signal(SIGFPE, catch_sigfpe);

  // Toggle the state.
  fpeTrappingActive = true;
  return fpeTrappingActive;
}

//------------------------------------------------------------------------------------------------//
//! Disable trapping of floating point errors.
void fpe_trap::disable() {
  ieee_set_fp_control(0x00);
  fpeTrappingActive = false;
  return;
}

} // end namespace rtt_dsxx

#endif // FPETRAP_OSF_ALPHA

//------------------------------------------------------------------------------------------------//
// WINDOWS X86
//------------------------------------------------------------------------------------------------//
#ifdef FPETRAP_WINDOWS_X86

#include <Windows.h> // defines STATUS_FLOAT_...
#include <float.h>   // defines _controlfp_s
#include <intrin.h>  // _ReturnAddress
#include <new.h>     // _set_new_handler
#include <signal.h>  // SIGABRT

// typdef ignored on left...
#pragma warning(push)
#pragma warning(disable : 4091)
#include <dbghelp.h> // minidump_exception_information
#pragma warning(pop)

#if !defined(__clang__)
#pragma fenv_access(on)
#endif

/* Signal handler for floating point exceptions. */
extern "C" void trans_func(unsigned int u, PEXCEPTION_POINTERS /*pExp*/) {
  std::cout << "(fpe_trap/windows_x86.cc) A SIGFPE was detected!" << std::endl;

  std::string mesg;
  switch (u) {
  case STATUS_INTEGER_DIVIDE_BY_ZERO:
    mesg = "Integer divide by zero.";
    break;
  case STATUS_INTEGER_OVERFLOW:
    mesg = "Integer overflow.";
    break;
  case STATUS_FLOAT_DIVIDE_BY_ZERO:
    mesg = "Floating point divide by zero.";
    break;
  case STATUS_FLOAT_OVERFLOW:
    mesg = "Floating point overflow.";
    break;
  case STATUS_FLOAT_UNDERFLOW:
    mesg = "Floating point underflow.";
    break;
  case STATUS_FLOAT_INEXACT_RESULT:
    mesg = "Floating point inexact result.";
    break;
  case STATUS_FLOAT_INVALID_OPERATION:
    mesg = "Invalid floating point operation.";
    break;
  default:
    mesg = "Unknown floating point exception.";
    break;
  }

  Insist(0, mesg);
}

namespace rtt_dsxx {

// ------------------------------------------------------------------------------------------------
// - http://stackoverflow.com/questions/2769814/how-do-i-use-try-catch-to-catch-floating-point-errors
// - See MSDN articles on fenv_access and _controlfp_s examples.
// ------------------------------------------------------------------------------------------------
bool fpe_trap::enable() {
  // Always call this before setting control words.
  _clearfp();

  // Read the current control words.
  unsigned int fp_control_word = _controlfp(0, 0);

  // Set the exception masks off for exceptions that you want to trap.  When a mask bit is set, the
  // corresponding floating-point exception is blocked from being generated.

  fp_control_word &= ~(EM_INVALID | EM_ZERODIVIDE | EM_OVERFLOW);

  // Other options: _EM_DENORMAL, _EM_UNDERFLOW, _EM_INEXACT

  // Update the control word with our changes: MCW_EM is Interrupt exception mask.
  _controlfp(fp_control_word, MCW_EM);

  if (this->abortWithInsist)
    _set_se_translator(trans_func);

  // Toggle the state.
  fpeTrappingActive = true;

  return fpeTrappingActive;
}

//------------------------------------------------------------------------------------------------//
//! Disable trapping of floating point errors.
void fpe_trap::disable() {
  // Always call this before setting control words.
  _clearfp();

  // Read the current control words.
  unsigned int fp_control_word = _controlfp(0, 0);

  // Set the exception masks off for exceptions that you want to trap.  When a mask bit is set, the
  // corresponding floating-point exception is blocked from being generated.
  fp_control_word |= (EM_INVALID | EM_ZERODIVIDE | EM_OVERFLOW);

  // Update the control word with our changes MCW_EM is Interrupt exception mask.
  _controlfp(fp_control_word, MCW_EM);

  fpeTrappingActive = false;
  return;
}

//------------------------------------------------------------------------------------------------//
// CCrashHandler
//------------------------------------------------------------------------------------------------//
void CCrashHandler::SetProcessExceptionHandlers() {
  std::cout << "In CCrashHandler::SetProcessExceptionHandlers" << std::endl;
  // Install top-level SEH handler
  SetUnhandledExceptionFilter(SehHandler);

  // Catch pure virtual function calls.  Because there is one _purecall_handler for the whole
  // process, calling this function immediately impacts all threads. The last caller on any thread
  // sets the handler.  http://msdn.microsoft.com/en-us/library/t296ys27.aspx
  _set_purecall_handler(PureCallHandler);

  // Catch new operator memory allocation exceptions
  _set_new_handler(NewHandler);

  // Catch invalid parameter exceptions.
  _set_invalid_parameter_handler(InvalidParameterHandler);

  // Set up C++ signal handlers

  _set_abort_behavior(_CALL_REPORTFAULT, _CALL_REPORTFAULT);

  // Catch an abnormal program termination
  signal(SIGABRT, SigabrtHandler);

  // Catch illegal instruction handler
  signal(SIGINT, SigintHandler);

  // Catch a termination request
  signal(SIGTERM, SigtermHandler);
}

void CCrashHandler::SetThreadExceptionHandlers() {
  // Catch terminate() calls.  In a multithreaded environment, terminate functions are maintained
  // separately for each thread. Each new thread needs to install its own terminate function. Thus,
  // each thread is in charge of its own termination handling.
  // http://msdn.microsoft.com/en-us/library/t6fk7h29.aspx
  set_terminate(TerminateHandler);

  // Catch unexpected() calls.  In a multithreaded environment, unexpected functions are maintained
  // separately for each thread. Each new thread needs to install its own unexpected function. Thus,
  // each thread is in charge of its own unexpected handling.
  // http://msdn.microsoft.com/en-us/library/h46t5b69.aspx
  set_unexpected(UnexpectedHandler);

  // Catch a floating point error
  typedef void (*sigh)(int);
  signal(SIGFPE, (sigh)SigfpeHandler);

  // Catch an illegal instruction
  signal(SIGILL, SigillHandler);

  // Catch illegal storage access errors
  signal(SIGSEGV, SigsegvHandler);
}

// The following code gets exception pointers using a workaround found in CRT code.
void CCrashHandler::GetExceptionPointers(DWORD dwExceptionCode,
                                         EXCEPTION_POINTERS **ppExceptionPointers) {
  // The following code was taken from VC++ 8.0 CRT (invarg.c: line 104)

  EXCEPTION_RECORD ExceptionRecord;
  CONTEXT ContextRecord;
  memset(&ContextRecord, 0, sizeof(CONTEXT));

#ifdef _X86_

  __asm {
        mov dword ptr [ContextRecord.Eax], eax
        mov dword ptr [ContextRecord.Ecx], ecx
        mov dword ptr [ContextRecord.Edx], edx
        mov dword ptr [ContextRecord.Ebx], ebx
        mov dword ptr [ContextRecord.Esi], esi
        mov dword ptr [ContextRecord.Edi], edi
        mov word ptr [ContextRecord.SegSs], ss
        mov word ptr [ContextRecord.SegCs], cs
        mov word ptr [ContextRecord.SegDs], ds
        mov word ptr [ContextRecord.SegEs], es
        mov word ptr [ContextRecord.SegFs], fs
        mov word ptr [ContextRecord.SegGs], gs
        pushfd
        pop [ContextRecord.EFlags]
  }

  ContextRecord.ContextFlags = CONTEXT_CONTROL;
#pragma warning(push)
#pragma warning(disable : 4311)
  ContextRecord.Eip = (ULONG)_ReturnAddress();
  ContextRecord.Esp = (ULONG)_AddressOfReturnAddress();
#pragma warning(pop)
  ContextRecord.Ebp = *((ULONG *)_AddressOfReturnAddress() - 1);

#elif defined(_IA64_) || defined(_AMD64_)

  /* Need to fill up the Context in IA64 and AMD64. */
  RtlCaptureContext(&ContextRecord);

#else /* defined (_IA64_) || defined (_AMD64_) */

  ZeroMemory(&ContextRecord, sizeof(ContextRecord));

#endif /* defined (_IA64_) || defined (_AMD64_) */

  ZeroMemory(&ExceptionRecord, sizeof(EXCEPTION_RECORD));

  ExceptionRecord.ExceptionCode = dwExceptionCode;
  ExceptionRecord.ExceptionAddress = _ReturnAddress();

  ///

  EXCEPTION_RECORD *pExceptionRecord = new EXCEPTION_RECORD;
  memcpy(pExceptionRecord, &ExceptionRecord, sizeof(EXCEPTION_RECORD));
  CONTEXT *pContextRecord = new CONTEXT;
  memcpy(pContextRecord, &ContextRecord, sizeof(CONTEXT));

  *ppExceptionPointers = new EXCEPTION_POINTERS;
  (*ppExceptionPointers)->ExceptionRecord = pExceptionRecord;
  (*ppExceptionPointers)->ContextRecord = pContextRecord;
}

//! Action to perform when an exception is found.
void CCrashHandler::ActionOnException(std::string const &message, EXCEPTION_POINTERS *pExcPtrs) {
  std::cout << "In CCrashHandler::ActionOnException" << std::endl;
  CCrashHandler::CreateMiniDump(pExcPtrs);

  // Provide a report
  std::ostringstream msg;
  msg << "\nERROR: " << message << "\nA fatal error has occured and a dump file was generated "
      << "(crashdump.dmp)."
      << "\nThe dump file can be viewed by loading it into "
      << "Visual Studio." << std::endl;
  std::cout << msg.str() << std::endl;

  // Terminate process
  TerminateProcess(GetCurrentProcess(), 1);

  return;
}

// This method creates minidump of the process
void CCrashHandler::CreateMiniDump(EXCEPTION_POINTERS *pExcPtrs) {
  std::cout << "In CCrashHandler::CreateMiniDump" << std::endl;
#define _T(foo) foo
  HMODULE hDbgHelp = NULL;
  HANDLE hFile = NULL;
  MINIDUMP_EXCEPTION_INFORMATION mei;
  MINIDUMP_CALLBACK_INFORMATION mci;

  // Load dbghelp.dll
  hDbgHelp = LoadLibrary(_T("dbghelp.dll"));
  if (hDbgHelp == NULL)
    return;

  // Create the minidump file
  hFile = CreateFile(_T("crashdump.dmp"), GENERIC_WRITE, 0, NULL, CREATE_ALWAYS,
                     FILE_ATTRIBUTE_NORMAL, NULL);

  if (hFile == INVALID_HANDLE_VALUE)
    return;

  // Write minidump to the file
  mei.ThreadId = GetCurrentThreadId();
  mei.ExceptionPointers = pExcPtrs;
  mei.ClientPointers = FALSE;
  mci.CallbackRoutine = NULL;
  mci.CallbackParam = NULL;

  std::cout << "\nThreadId = " << mei.ThreadId << "\nExceptionPointers = " << mei.ExceptionPointers;

  typedef BOOL(WINAPI * LPMINIDUMPWRITEDUMP)(
      HANDLE hProcess, DWORD ProcessId, HANDLE hFile, MINIDUMP_TYPE DumpType,
      CONST PMINIDUMP_EXCEPTION_INFORMATION ExceptionParam,
      CONST PMINIDUMP_USER_STREAM_INFORMATION UserEncoderParam,
      CONST PMINIDUMP_CALLBACK_INFORMATION CallbackParam);

  LPMINIDUMPWRITEDUMP pfnMiniDumpWriteDump =
      (LPMINIDUMPWRITEDUMP)GetProcAddress(hDbgHelp, "MiniDumpWriteDump");
  if (!pfnMiniDumpWriteDump)
    return;

  HANDLE hProcess = GetCurrentProcess();
  DWORD dwProcessId = GetCurrentProcessId();

  std::cout << "\nhProcess = " << hProcess << "\nProcessId = " << dwProcessId
            << "\nhFile = " << hFile << "\nMiniDumpNormal = " << MiniDumpNormal;

  BOOL bWriteDump =
      pfnMiniDumpWriteDump(hProcess, dwProcessId, hFile, MiniDumpNormal, &mei, NULL, &mci);

  if (!bWriteDump)
    return;

  // Close file
  CloseHandle(hFile);

  // Unload dbghelp.dll
  FreeLibrary(hDbgHelp);
}

// Structured exception handler
LONG WINAPI CCrashHandler::SehHandler(PEXCEPTION_POINTERS pExceptionPtrs) {
  std::cout << "In CCrashHandler::SehHandler" << std::endl;

  // Write minidump file
  ActionOnException(std::string("Seh"), pExceptionPtrs);

  // Unreacheable code
  return EXCEPTION_EXECUTE_HANDLER;
}

// CRT terminate() call handler
void __cdecl CCrashHandler::TerminateHandler() {
  std::cout << "In CCrashHandler::TerminateHandler" << std::endl;
  // Abnormal program termination (terminate() function was called)

  // Retrieve exception information
  EXCEPTION_POINTERS *pExceptionPtrs = NULL;
  GetExceptionPointers(0, &pExceptionPtrs);

  // Write minidump file
  ActionOnException(std::string("Terminate"), pExceptionPtrs);

  return;
}

// CRT unexpected() call handler
void __cdecl CCrashHandler::UnexpectedHandler() {
  std::cout << "In CCrashHandler::UnexpectedHandler" << std::endl;
  // Unexpected error (unexpected() function was called)

  // Retrieve exception information
  EXCEPTION_POINTERS *pExceptionPtrs = NULL;
  GetExceptionPointers(0, &pExceptionPtrs);

  // Write minidump file
  ActionOnException(std::string("Unexpected"), pExceptionPtrs);

  return;
}

// CRT Pure virtual method call handler
void __cdecl CCrashHandler::PureCallHandler() {
  std::cout << "In CCrashHandler::PurceCallHandler" << std::endl;
  // Pure virtual function call

  // Retrieve exception information
  EXCEPTION_POINTERS *pExceptionPtrs = NULL;
  GetExceptionPointers(0, &pExceptionPtrs);

  // Write minidump file
  ActionOnException(std::string("PureCall"), pExceptionPtrs);

  return;
}

// CRT invalid parameter handler
void __cdecl CCrashHandler::InvalidParameterHandler(const wchar_t * /*expression*/,
                                                    const wchar_t * /*function*/,
                                                    const wchar_t * /*file*/, unsigned int /*line*/,
                                                    uintptr_t /*pReserved*/) {

  std::cout << "In CCrashHandler::InvalidParameterHandler" << std::endl;
  // pReserved;

  // Invalid parameter exception

  // Retrieve exception information
  EXCEPTION_POINTERS *pExceptionPtrs = NULL;
  GetExceptionPointers(0, &pExceptionPtrs);

  // Write minidump file
  ActionOnException(std::string("InvalidParameter"), pExceptionPtrs);

  return;
}

// CRT new operator fault handler
int __cdecl CCrashHandler::NewHandler(size_t) {
  std::cout << "In CCrashHandler::NewHandler" << std::endl;
  // 'new' operator memory allocation exception

  // Retrieve exception information
  EXCEPTION_POINTERS *pExceptionPtrs = NULL;
  GetExceptionPointers(0, &pExceptionPtrs);

  // Write minidump file
  ActionOnException(std::string("MemoryAllocation"), pExceptionPtrs);

  return 0;
}

// CRT SIGABRT signal handler
void CCrashHandler::SigabrtHandler(int) {
  std::cout << "In CCrashHandler::SigabrtHandler" << std::endl;
  // Caught SIGABRT C++ signal

  // Retrieve exception information
  EXCEPTION_POINTERS *pExceptionPtrs = NULL;
  GetExceptionPointers(0, &pExceptionPtrs);

  // Write minidump file
  ActionOnException(std::string("SignalAbort"), pExceptionPtrs);

  return;
}

// CRT SIGFPE signal handler
void CCrashHandler::SigfpeHandler(int /*code*/, int /*subcode*/) {
  std::cout << "In CCrashHandler::SigfpeHandler" << std::endl;
  // Floating point exception (SIGFPE)

  EXCEPTION_POINTERS *pExceptionPtrs = (PEXCEPTION_POINTERS)_pxcptinfoptrs;

  // Write minidump file
  ActionOnException(std::string("SignalFloatingPointException"), pExceptionPtrs);

  return;
}

// CRT sigill signal handler
void CCrashHandler::SigillHandler(int) {
  std::cout << "In CCrashHandler::SigillHandler" << std::endl;
  // Illegal instruction (SIGILL)

  // Retrieve exception information
  EXCEPTION_POINTERS *pExceptionPtrs = NULL;
  GetExceptionPointers(0, &pExceptionPtrs);

  // Write minidump file
  ActionOnException(std::string("SignalIllegalInstruction"), pExceptionPtrs);

  return;
}

// CRT sigint signal handler
void CCrashHandler::SigintHandler(int) {
  std::cout << "In CCrashHandler::SigintHandler" << std::endl;
  // Interruption (SIGINT)

  // Retrieve exception information
  EXCEPTION_POINTERS *pExceptionPtrs = NULL;
  GetExceptionPointers(0, &pExceptionPtrs);

  // Write minidump file
  ActionOnException(std::string("SignalInterruption"), pExceptionPtrs);

  return;
}

// CRT SIGSEGV signal handler
void CCrashHandler::SigsegvHandler(int) {
  std::cout << "In CCrashHandler::SigsegvHandler" << std::endl;
  // Invalid storage access (SIGSEGV)

  PEXCEPTION_POINTERS pExceptionPtrs = (PEXCEPTION_POINTERS)_pxcptinfoptrs;

  // Write minidump file
  ActionOnException(std::string("SignalInvalidStorageAccess(SIGSEGV)"), pExceptionPtrs);

  return;
}

// CRT SIGTERM signal handler
void CCrashHandler::SigtermHandler(int) {
  std::cout << "In CCrashHandler::SigtermHandler" << std::endl;
  // Termination request (SIGTERM)

  // Retrieve exception information
  EXCEPTION_POINTERS *pExceptionPtrs = NULL;
  GetExceptionPointers(0, &pExceptionPtrs);

  // Write minidump file
  ActionOnException(std::string("SignalTerminateRequest"), pExceptionPtrs);

  return;
}

} // end namespace rtt_dsxx

#endif // FPETRAP_WINDOWS_X86

//------------------------------------------------------------------------------------------------//
// DARWIN INTEL
//------------------------------------------------------------------------------------------------//
#ifdef FPETRAP_DARWIN_INTEL

#include <xmmintrin.h>

namespace rtt_dsxx {

//------------------------------------------------------------------------------------------------//
//! Enable trapping of floating point errors.
bool fpe_trap::enable() {
  _mm_setcsr(_MM_MASK_MASK & ~(_MM_MASK_OVERFLOW | _MM_MASK_INVALID | _MM_MASK_DIV_ZERO));

  // This functionality is not currently implemented on the Mac.
  abortWithInsist = false;

  // Toggle the state.
  fpeTrappingActive = true;
  return fpeTrappingActive;
}
//------------------------------------------------------------------------------------------------//
//! Disable trapping of floating point errors.
void fpe_trap::disable() {
  _mm_setcsr(_MM_MASK_MASK & ~0x00);
  fpeTrappingActive = false;
  return;
}

} // namespace rtt_dsxx

#endif // FPETRAP_DARWIN_INTEL

//------------------------------------------------------------------------------------------------//
// DARWIN PPC
//------------------------------------------------------------------------------------------------//
#ifdef FPETRAP_DARWIN_PPC

#include <mach/mach.h>

// Local functions

namespace {

/*
 * On Mach, we need a mindbogglingly complex setup for floating point errors.  Not the least of the
 * hassles is that we have to do the whole thing from a different thread.
 */
void *fpe_enabler(void *parent) {
  mach_port_t victim = (mach_port_t)parent;
  mach_msg_type_number_t count;

  ppc_thread_state_t ts;
  ppc_float_state_t fs;

  /* First enable the right FP exception conditions */
  count = PPC_FLOAT_STATE_COUNT;
  thread_get_state(victim, PPC_FLOAT_STATE, (thread_state_t)&fs, &count);
  /* Enable VE OE ZE, Disable UE XE */
  fs.fpscr = (fs.fpscr & ~0x1FFFFF28) | 0x0D0;
  thread_set_state(victim, PPC_FLOAT_STATE, (thread_state_t)&fs, count);

  /* Now, enable FP exceptions as such */
  count = PPC_THREAD_STATE_COUNT;
  thread_get_state(victim, PPC_THREAD_STATE, (thread_state_t)&ts, &count);
  /* Set FE0 = FE1 = 1 */
  ts.srr1 |= 0x0900;
  thread_set_state(victim, PPC_THREAD_STATE, (thread_state_t)&ts, count);

  return 0;
}

static void catch_sigfpe(int sig) { Insist(0, "Floating point exception caught by fpe_trap.") }

} // end of namespace

namespace rtt_dsxx {
//!  Enable trapping of floating point errors.
bool fpe_trap::enable() {
  pthread_t enabler;
  void *mts = reinterpret_cast<void *>(mach_thread_self());
  pthread_create(&enabler, NULL, fpe_enabler, mts);
  pthread_join(enabler, NULL);

  if (this->abortWithInsist)
    signal(SIGFPE, catch_sigfpe);

  // Toggle the state.
  fpeTrappingActive = true;
  return fpeTrappingActive;
}

//------------------------------------------------------------------------------------------------//
//! Disable trapping of floating point errors.
void fpe_trap::disable() {
  // (void)feenableexcept( 0x00 );
  Insist(0, "Please update darwin_ppc.cc to provide instructions for disabling "
            "fpe traps.");
  // fpeTrappingActive=false;
  return;
}

} // namespace rtt_dsxx

#endif // FPETRAP_DARWIN_PPC

//------------------------------------------------------------------------------------------------//
// UNSUPPORTED PLATFORMS
//------------------------------------------------------------------------------------------------//
#ifdef FPETRAP_UNSUPPORTED

namespace rtt_dsxx {

//------------------------------------------------------------------------------------------------//
//!  Enable trapping of floating point errors.
bool fpe_trap::enable() {
  // unsupported platform.  leave fag set to false.  (using abortWithInsist to silence unused
  // variable warning)
  if (abortWithInsist)
    fpeTrappingActive = false;
  else
    fpeTrappingActive = false;

  return fpeTrappingActive;
}

//------------------------------------------------------------------------------------------------//
//! Disable trapping of floating point errors.
void fpe_trap::disable() { return; }

} // end namespace rtt_dsxx

#endif // FPETRAP_UNSUPPORTED

//------------------------------------------------------------------------------------------------//
// end of ds++/fpe_trap.cc
//------------------------------------------------------------------------------------------------//
