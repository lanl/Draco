//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   diagnostics/Diagnostics.hh
 * \author Thomas M. Evans, Aimee Hungerford
 * \date   Fri Dec  9 10:52:38 2005
 * \brief  Diagnostics class for runtime info.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef diagnostics_Diagnostics_hh
#define diagnostics_Diagnostics_hh

#include "diagnostics/config.h"
#include "ds++/config.h"
#include <map>
#include <string>
#include <vector>

namespace rtt_diagnostics {

//================================================================================================//
/*!
 * \namespace Diagnostics
 * \brief Allows code clients to register diagnostics data during runtime.
 *
 * This namespace defines maps for the following types:
 * - int
 * - double
 * - vector<int>
 * - vector<double>
 * .
 * The key for each map is a std::string.
 *
 * These maps can be used to store diagnostic quantities.  Because they have global, file scope they
 * can be accessed from any routine.  The general usage is as follows:
 *
 * \code
 *   rtt_diagnostics::Diagnostics::integers["Num_part_per_proc"] = 1011;
 * \endcode
 *
 * A compile-time switching mechanism for using these maps is provided by the macros
 * \c DIAGNOSTICS_ONE, \c DIAGNOSTICS_TWO, and \c DIAGNOSTICS_THREE.
 *
 * \example diagnostics/test/tstDiagnostics.cc
 */
//================================================================================================//

namespace Diagnostics {

//! Map of integer data.
extern DLL_PUBLIC_diagnostics std::map<std::string, int> integers;

//! Map of floating point data.
extern DLL_PUBLIC_diagnostics std::map<std::string, double> doubles;

//! Map of vector, integer data.
extern DLL_PUBLIC_diagnostics std::map<std::string, std::vector<int>> vec_integers;

//! Map of vector, double data.
extern DLL_PUBLIC_diagnostics std::map<std::string, std::vector<double>> vec_doubles;

} // end of namespace Diagnostics

} // end namespace rtt_diagnostics

//------------------------------------------------------------------------------------------------//
/*!
 * \page diagnostics Diagnostics Levels
 *
 * The diagnostics can be turned on in three different levels based on logical bit comparisons.
 * The following shows the levels: 
 * - Bit 0, (001), activates Level 1 (negligible performance impact)
 * - Bit 1, (010), activates Level 2 (some performance impact and possible intrusive output, 
 *                 rtt_memory tracking is activated.)
 * - Bit 2, (100), activates Level 3 (includes fpe_trap diagnostics)
 * .
 * The following integer settings activate Levels in the following way:
 * - 0 all off
 * - 1 Level 1
 * - 2 Level 2
 * - 3 Level 1, Level 2
 * - 4 Level 3
 * - 5 Level 1, Level 3
 * - 6 Level 2, Level 3
 * - 7 Level 1, Level 2, Level 3
 * .
 * Thus setting \c -DDRACO_DIAGNOSTICS=7 at configure time will turn on all levels.  The default
 * setting is 0.
 *
 * The intent is to use Level 1 for high-level, low cost diagnostics that are always active (ie.
 * User "Education").  Levels 2 and 3 are for low-level diagnostics that could incur a performance
 * penalty.  However, all of these usages are up to the client.
 *
 * The value for \c DRACO_DIAGNOSTICS is set and saved in ds++'s \c CMakeLists.txt and \c config.h,
 * respectively.
 */
/*!
 * \def DIAGNOSTICS_ONE(Diagnostics::member)
 *
 * Single-line statement macro for diagnostics Level 1:
 * \code
 *     DIAGNOSTICS_ONE(integers["Variable"] = 1);
 * \endcode
 * On when DRACO_DIAGNOSTICS & 1 is true.  Defines DRACO_DIAGNOSTICS_LEVEL_1.
 */
/*!
 * \def DIAGNOSTICS_TWO(Diagnostics::member)
 *
 * Single-line statement macro for diagnostics Level 2:
 * \code
 *     DIAGNOSTICS_TWO(integers["Variable"] = 1);
 * \endcode
 * On when DRACO_DIAGNOSTICS & 2 is true.  Defines DRACO_DIAGNOSTICS_LEVEL_2.
 */
/*!
 * \def DIAGNOSTICS_THREE(Diagnostics::member)
 *
 * Single-line statement macro for diagnostics Level 3:
 * \code
 *     DIAGNOSTICS_THREE(integers["Variable"] = 1);
 * \endcode
 * On when DRACO_DIAGNOSTICS & 4 is true.  Defines DRACO_DIAGNOSTICS_LEVEL_3.
 */
//------------------------------------------------------------------------------------------------//
#ifdef DRACO_DIAGNOSTICS_LEVEL_1
#define DIAGNOSTICS_ONE(member) rtt_diagnostics::Diagnostics::member
#else
#define DIAGNOSTICS_ONE(member)
#endif

#ifdef DRACO_DIAGNOSTICS_LEVEL_2
#define DIAGNOSTICS_TWO(member) rtt_diagnostics::Diagnostics::member
#else
#define DIAGNOSTICS_TWO(member)
#endif

#ifdef DRACO_DIAGNOSTICS_LEVEL_3
#define DIAGNOSTICS_THREE(member) rtt_diagnostics::Diagnostics::member
#else
#define DIAGNOSTICS_THREE(member)
#endif

#endif // diagnostics_Diagnostics_hh

//------------------------------------------------------------------------------------------------//
// end of diagnostics/Diagnostics.hh
//------------------------------------------------------------------------------------------------//
