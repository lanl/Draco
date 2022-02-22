//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   diagnostics/Procmon.hh
 * \author Kelly Thompson
 * \date   Monday, Apr 22, 2013, 10:10 am
 * \brief  Procmon class for printing runtime system diagnostics (free memory per node, etc).
 * \note   Copyright (C) 2013-2022 Triad National Security, LLC., All rights reserved.
 *
 * A file with the same name lives in wedgehog_components and provided the macros \c PROCMON_START
 * and \c PROCMON_STOP.  These command would provide a one line system summary (including timing
 * information).  We don't need the timing information as that can be accessed using other tools
 * (\em Draco timers or system profilers like \em HPC \em Toolkit).  This file provides a similar
 * one line diagnostic output that provides a snapshot in time for system memory via the Linux
 * system command \c 'getrusage.'
 *
 * Mike McKay provides a package in \em xRage that provides these memory instrumentation calls and
 * more for many systems.  It may be useful to review his coding when extending the feature set of
 * this package or when porting these capabilities to new architectures (i.e.: IBM BlueGene ).  His
 * files are: \c procmon_info.h, \c procmon_info.c, \c procmon_iface.h and \c resources.f90. */
//------------------------------------------------------------------------------------------------//

#ifndef diagnostics_Procmon_hh
#define diagnostics_Procmon_hh

#include "diagnostics/config.h"
#include "ds++/config.h"
#include <iostream>
#include <string>
#include <vector>

namespace rtt_diagnostics {

//================================================================================================//
/*!
 * \brief Access and report on system and process resources.
 *
 * \param[in] identifier
 * \param[in] mynode
 * \return    void
 *
 * Examine /proc/pid/status and /proc/meminfo and report results.
 *
 * Sample meminfo content:
 *
 *     MemTotal:     32949476 kB
 *     MemFree:      25156352 kB
 *     Buffers:        525952 kB
 *     Cached:        5777456 kB
 *     SwapCached:      65124 kB
 *     Active:        3844364 kB
 *     Inactive:      2853604 kB
 *     HighTotal:           0 kB
 *     HighFree:            0 kB
 *     LowTotal:     32949476 kB
 *     LowFree:      25156352 kB
 *     SwapTotal:    34996216 kB
 *     SwapFree:     34884456 kB
 *     Dirty:              16 kB
 *     Writeback:           0 kB
 *     AnonPages:      382284 kB
 *     Mapped:          33868 kB
 *     Slab:          1029048 kB
 *     PageTables:      13592 kB
 *     NFS_Unstable:        0 kB
 *     Bounce:              0 kB
 *     CommitLimit:  51470952 kB
 *     Committed_AS:   812588 kB
 *     VmallocTotal: 34359738367 kB
 *     VmallocUsed:    280384 kB
 *     VmallocChunk: 34359456351 kB
 *     HugePages_Total:     0
 *     HugePages_Free:      0
 *     HugePages_Rsvd:      0
 *     Hugepagesize:     2048 kB
 *
 * Sample /proc/pid/status content:
 *
 *     Name:   tstProcmon
 *     State:  R (running)
 *     SleepAVG:       98%
 *     Tgid:   14768
 *     Pid:    14768
 *     PPid:   14767
 *     TracerPid:      0
 *     Uid:    2017    2017    2017    2017
 *     Gid:    2017    2017    2017    2017
 *     FDSize: 64
 *     Groups: 217 1462 2017 5016 8357 17872 18662 19038 21356 22019 22071
 *     VmPeak:   106392 kB
 *     VmSize:   102084 kB
 *     VmLck:         0 kB
 *     VmHWM:      4616 kB
 *     VmRSS:      4616 kB
 *     VmData:     3000 kB
 *     VmStk:        92 kB
 *     VmExe:        16 kB
 *     VmLib:      6316 kB
 *     VmPTE:       268 kB
 *     StaBrk: 0bdb0000 kB
 *     Brk:    0be5a000 kB
 *     StaStk: 7fff0a0c22b0 kB
 *     Threads:        1
 *     SigQ:   1/270336
 *     SigPnd: 0000000000000000
 *     ShdPnd: 0000000000000000
 *     SigBlk: 0000000000000000
 *     SigIgn: 0000000000000000
 *     SigCgt: 00000001800104e0
 *     CapInh: 0000000000000000
 *     CapPrm: 0000000000000000
 *     CapEff: 0000000000000000
 *     Cpus_allowed:   00000000,00000000,00000000,00000000,00000000,00000000,00000000,ffffffff
 *     Mems_allowed:   00000000,00000003
 *
 * \example diagnostics/test/tstProcmon.cc
 */
//================================================================================================//
void procmon_resource_print(std::string const &identifier, int const &mynode = -1,
                            std::ostream &msg = std::cout);

} // end namespace rtt_diagnostics

#ifdef USE_PROCMON
#define PROCMON_REPORT(string) rtt_diagnostics::procmon_resource_print(string)
#else
#define PROCMON_REPORT(string)
#endif

#endif // diagnostics_Procmon_hh

//------------------------------------------------------------------------------------------------//
// end of diagnostics/Procmon.hh
//------------------------------------------------------------------------------------------------//
