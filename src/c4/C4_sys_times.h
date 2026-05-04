//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   c4/C4_sys_times.h
 * \author Kelly Thompson
 * \brief  Encapsulate system headers for timing information.
 * \note   Copyright (C) 2010-2026 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef c4_C4_sys_times_h
#define c4_C4_sys_times_h

#include "dsxx/config.h"

#if defined(_MSC_VER) || defined(MINGW)
#include <chrono>
#include <ctime>
#else
#include <sys/times.h>
#include <unistd.h>
#endif

#if defined(_MSC_VER)
// Consider using GetProcessTimes instead of this as a higher resolution timer.
#define DRACO_TIME_TYPE std::chrono::high_resolution_clock::time_point
#define DRACO_CLOCKS_PER_SEC CLOCKS_PER_SEC
#else

#define DRACO_TIME_TYPE tms
#define DRACO_CLOCKS_PER_SEC sysconf(_SC_CLK_TCK)
#endif

#endif // c4_C4_sys_times_h

//------------------------------------------------------------------------------------------------//
// end of c4/C4_sys_times.h
//------------------------------------------------------------------------------------------------//
