/*---------------------------------------------*-C-*----------------------------------------------*/
/*!
 * \file   parser/Debug_Options.hh
 * \author Kent Grimmett Budge
 * \brief  Define the Debug_Options enumeration and declare parse functions.
 * \note   Copyright (C) 2016-2022 Triad National Security, LLC., All rights reserved. */
/*------------------------------------------------------------------------------------------------*/

#ifndef parser_Debug_Options_hh
#define parser_Debug_Options_hh

#include "Token_Stream.hh"

namespace rtt_parser {
//------------------------------------------------------------------------------------------------//
//! Enumeration of debug flag bits

enum Debug_Options : unsigned {
  DEBUG_NONE = 0,           // For human readability; No debug options desired.
  DEBUG_SUMMARY = 1,        // Report *short* summary of balance/performance.
  DEBUG_BALANCE = 2,        // Report on energy balance.
  DEBUG_TIMESTEP = 4,       // Report on what is limiting the time step.
  DEBUG_TIMING = 8,         // Report CPU times for various code regions.
  DEBUG_ALGORITHM = 16,     // Report on behavior of the algorithm.
  DEBUG_MEMORY = 32,        // Report on memory usage.
  DEBUG_PROBLEM = 64,       // Report on the characteristics of the posed problem; e.g. what is its
                            // scattering fraction?
  DEBUG_GMV_DUMP = 128,     // Produce a GMV dump of the solution.
  DEBUG_RESET_TIMING = 256, // Reset all timings to zero.

  DEBUG_END = 512, // Sentinel value and first available extension.

  DEBUG_SILENT_MASK = ~(DEBUG_RESET_TIMING)
  // Options producing no terminal output
};

//------------------------------------------------------------------------------------------------//
//! Parse debug options in uniform way
unsigned parse_debug_options(rtt_parser::Token_Stream &tokens, unsigned parent_mask = 0);

//------------------------------------------------------------------------------------------------//
//! Add an application-specific debug option.
void add_debug_option(string const &option_name, unsigned const bit);

//------------------------------------------------------------------------------------------------//
//! Add an application-specific debug option.
unsigned add_debug_option(string const &option_name);

//------------------------------------------------------------------------------------------------//
unsigned get_debug_option(string const &option_name);

//------------------------------------------------------------------------------------------------//
//! Flush application-specific debug options.
void flush_debug_options();

//------------------------------------------------------------------------------------------------//
//! Write debug options in a manner that can be parsed
std::string debug_options_as_text(unsigned debug_options);

} // namespace rtt_parser

#endif /* parser_Debug_Options_hh */

/*------------------------------------------------------------------------------------------------*/
/* end of parser/Debug_Options.hh */
/*------------------------------------------------------------------------------------------------*/
