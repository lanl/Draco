/*---------------------------------------------*-C-*----------------------------------------------*/
/*!
 * \file   parser/Debug_Options.cc
 * \author Kent Grimmett Budge
 * \brief  Define Debug_Options parse functions.
 * \note   Copyright (C) 2016-2022 Triad National Security, LLC., All rights reserved. */
/*------------------------------------------------------------------------------------------------*/

#include "Debug_Options.hh"
#include <map>

namespace // anonymous
{

unsigned available = rtt_parser::DEBUG_END;

std::map<std::string, unsigned> extended_debug_option;
std::map<unsigned, std::string> extended_debug_back_option;

#ifdef REQUIRE_ON

//! Is the bit actually a bit? That is, a power of 2?
bool is_bit(unsigned bit) {
  Require(bit > 0); // corner case won't work

  // Shift to first bit
  while ((bit & 1U) == 0U) {
    bit >>= 1U;
  }
  // Erase that bit; see if the result is zero, as it must be for a power of 2.
  return (bit ^ 1U) == 0U;
}

#endif

} // end anonymous namespace

namespace rtt_parser {
using std::string;

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Get a debug specification.
 *
 * \param[in] option_name A string specifying a debug option keyword.
 * \return The bitmask value assigned to the keyword, or 0 if the keyword is not recognized.
 */
unsigned get_debug_option(string const &option_name) {
  if (option_name == "SUMMARY") {
    return DEBUG_SUMMARY;
  } else if (option_name == "BALANCE") {
    return DEBUG_BALANCE;
  } else if (option_name == "TIMESTEP") {
    return DEBUG_TIMESTEP;
  } else if (option_name == "TIMING") {
    return DEBUG_TIMING;
  } else if (option_name == "ALGORITHM") {
    return DEBUG_ALGORITHM;
  } else if (option_name == "MEMORY") {
    return DEBUG_MEMORY;
  } else if (option_name == "PROBLEM") {
    return DEBUG_PROBLEM;
  } else if (option_name == "GMV_DUMP") {
    return DEBUG_GMV_DUMP;
  } else if (option_name == "RESET_TIMING") {
    return DEBUG_RESET_TIMING;
  } else {
    // parse extension to debug options
    if (extended_debug_option.find(option_name) == extended_debug_option.end()) {
      return 0;
    } else {
      return extended_debug_option[option_name];
    }
  }
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Parse a debug specification.
 *
 * \param[in,out] tokens Token stream from which to parse a debug specification. The specification
 *      is a set of debug keywords, each optionally prefixed with a '!', and ends with the first
 *      token that is not a recognized debug keyword.
 * \param[in] parent Optional parent mask; defaults to zero. Allows a debug mask to be based on a
 *      parent mask, with selected bits added or masked out.
 * \return A debug options mask.
 */
unsigned parse_debug_options(Token_Stream &tokens, unsigned const parent) {
  using namespace rtt_parser;

  unsigned Result = parent;

  Token token = tokens.lookahead();
  while (token.type() == KEYWORD || token.text() == "!") {
    bool mask_out = (token.text() == "!");
    if (mask_out) {
      tokens.shift();
      token = tokens.lookahead();
    }
    unsigned mask = get_debug_option(token.text());
    if (mask != 0) {
      if (mask_out) {
        Result = Result & ~mask;
      } else {
        Result = Result | mask;
      }
      tokens.shift();
    } else {
      tokens.check_syntax(!mask_out, "trailing '!'");
      return Result;
    }
    token = tokens.lookahead();
  }
  return Result;
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Convert a debug mask to a string containing comma-delimited set of debug keywords.
 *
 * \param[in] debug_options Debug mask to be converted to a set of keywords.
 * \return A string containing a comma-delimited set of debug options.
 */
string debug_options_as_text(unsigned debug_options) {
  string Result;

  if (debug_options & DEBUG_SUMMARY) {
    Result += ", SUMMARY";
  }
  if (debug_options & DEBUG_BALANCE) {
    Result += ", BALANCE";
  }
  if (debug_options & DEBUG_TIMESTEP) {
    Result += ", TIMESTEP";
  }
  if (debug_options & DEBUG_TIMING) {
    Result += ", TIMING";
  }
  if (debug_options & DEBUG_ALGORITHM) {
    Result += ", ALGORITHM";
  }
  if (debug_options & DEBUG_MEMORY) {
    Result += ", MEMORY";
  }
  if (debug_options & DEBUG_PROBLEM) {
    Result += ", PROBLEM";
  }
  if (debug_options & DEBUG_GMV_DUMP) {
    Result += ", GMV_DUMP";
  }
  if (debug_options & DEBUG_RESET_TIMING) {
    Result += ", RESET_TIMING";
  }
  // Mask out standard options and see if any extensions are active
  debug_options =
      debug_options & ~(DEBUG_SUMMARY | DEBUG_BALANCE | DEBUG_TIMESTEP | DEBUG_TIMING |
                        DEBUG_ALGORITHM | DEBUG_MEMORY | DEBUG_PROBLEM | DEBUG_RESET_TIMING);

  if (debug_options) {
    for (const auto &i : extended_debug_option) {
      if (debug_options & i.second) {
        Result += ", " + i.first;
      }
    }
  }

  return Result;
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Add a new debug option to the debug parser specific to an application. This version
 *      assigns the next available bit.
 *
 * \param[in] option_name Debug option keyword
 * \return Bitflag value assigned to the new debug option.
 */
unsigned add_debug_option(string const &option_name) {
  if (extended_debug_option.find(option_name) != extended_debug_option.end()) {
    // option already exists; regard as benign
    return extended_debug_option[option_name];
  } else {
    while (available != 0 &&
           extended_debug_back_option.find(available) != extended_debug_back_option.end()) {
      available <<= 1U;
    }
    if (available == 0) {
      throw std::range_error("maximum debug options exceeded");
      // yeah, i know, if there are 4G debug options, someone has lost his
      // mind. Still.
    }
    extended_debug_option[option_name] = available;
    extended_debug_back_option[available] = option_name;
    return available;
  }
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Add a new debug option to the debug parser specific to an application. This version
 *      requests a specific bit and throws an exception if has already been requested elsewhere.
 *      This version will typically be called at the initial setup of an application.
 *
 * \param[in] option_name option keyword
 * \param[in] bit value to be assigned to the new debug option.
 */
void add_debug_option(string const &option_name, unsigned const bit) {
  Require(bit != 0);         // corner case will fail
  Require(bit >= DEBUG_END); // can't redefine standard debug
  Require(is_bit(bit));

  if (extended_debug_option.find(option_name) != extended_debug_option.end()) {
    if (extended_debug_option[option_name] != bit) {
      throw std::invalid_argument("debug option redefined");
    }
    // else duplicate identical definition acceptable
  } else if (extended_debug_back_option.find(bit) != extended_debug_back_option.end()) {
    throw std::invalid_argument("bitflag already allocated");
  } else {
    extended_debug_option[option_name] = bit;
    extended_debug_back_option[bit] = option_name;
  }
}

//------------------------------------------------------------------------------------------------//
void flush_debug_options() {
  extended_debug_option.clear();
  extended_debug_back_option.clear();
  available = DEBUG_END;
}

} // end namespace rtt_parser

/*------------------------------------------------------------------------------------------------*/
/* end of parser/Debug_Options.cc */
/*------------------------------------------------------------------------------------------------*/
